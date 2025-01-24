import csv
import json
import os
import re
import subprocess


ANNOTATION_PATH = "annotation.vcf.gz"
BGZIPPED_PATH = "bgzipped.bcf.gz"
HEADER_PATH = "header.vcf"
MAX_LINES_PER_PRINT = 100
MASK = "XXXXXXXXXX"

INFO_RESERVED_KEYS = {
    "AA",
    "AC",
    "AD",
    "ADF",
    "ADR",
    "AF",
    "AN",
    "BQ",
    "CIGAR",
    "DB",
    "DP",
    "END",
    "H2",
    "H3",
    "MQ",
    "MQ0",
    "NS",
    "SB",
    "SOMATIC",
    "VALIDATED",
    "1000G",
}

META_UNSTRUCTURED_WHITELIST = {
    "fileformat",
    "fileDate",
    "source",
    "reference",
    "assembly",
}

# Only check the Description value for these
META_STRUCTURED_WHITELIST = {
    "INFO",
    "FILTER",
    "FORMAT",
    "ALT",
    "contig",
}

PII_PATTERNS = [
    r"\b[a-zA-Z0-9._%+-]{3,}@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}\b",  # Email
    r"(?<!\.)\(?(?:\+62\)?[- ]?|\b0)(?:\d\)?[- ]?){1,3}\d{3,4}[- ]?\d{3,4}\b",  # Phone number
    r"(?<!SNOMED:)\b\d{16}\b",  # NIK
    r"\b(?:[A-Z]{1,2})\s?\d{1,4}\s?[A-Z]{0,3}\b",  # License plate
    r"\b(?:(?:Jl\.|Jalan|Desa|Kelurahan|Kecamatan|Kabupaten|Provinsi|Jakarta|Kode\s?Pos)(?:\s?(?:\d{5}|RT\s?\d{1,2}/RW\s?\d{1,2}|[A-Z^RT]+[a-z]*(?:\.\s?\d+)?),?)+,?\s?)+\b",  # Address
    r"\b(?:Dr\.|Prof\.|Ir\.|Haji|Hajjah|Putra|Putri|Sri|Adi)(?:\s[A-Z][a-z]+){1,2}\b",  # Name
]
ANY_PII_PATTERN = re.compile("|".join(f"(?:{pattern})" for pattern in PII_PATTERNS))

GENOMIC_SUFFIX_TYPES = {
    ".bcf": "u",
    ".bcf.gz": "b",
    ".vcf": "v",
    ".vcf.gz": "z",
}

METADATA_SUFFIXES = [
    ".json",
    ".csv",
    ".tsv",
    ".txt",
]

# Check if the script is running in AWS Lambda.
# EC2 instances don't have as much space in tmp
WORKING_DIR = os.getcwd()


class BcftoolsError(Exception):
    def __init__(self, message, stdout, stderr, returncode, process_args):
        self.message = message
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode
        self.process_args = process_args
        super().__init__(message)

    def __str__(self):
        return f"{self.message}\nProcess args: {self.process_args}\nstderr:\n{self.stderr}\nreturncode: {self.returncode}"


class ParsingError(Exception):
    pass


class CheckedProcess:
    def __init__(self, error_message, **kwargs):
        defaults = {
            "stderr": subprocess.PIPE,
            "cwd": WORKING_DIR,
            "encoding": "utf-8",
        }
        kwargs.update({k: v for k, v in defaults.items() if k not in kwargs})
        print(
            f"Running subprocess.Popen with kwargs: {json.dumps(kwargs, default=str)}"
        )
        self.process = subprocess.Popen(**kwargs)
        self.error_message = error_message
        self.stdout = self.process.stdout
        self.stdin = self.process.stdin

    def check(self):
        stdout, stderr = self.process.communicate()
        returncode = self.process.returncode
        if returncode != 0:
            raise BcftoolsError(
                self.error_message, stdout, stderr, returncode, self.process.args
            )


class Viewer:
    def __init__(self, file_path):
        self.started = False
        self.file_path = file_path
        self.view_process = None

    def _start(self):
        self.view_process = CheckedProcess(
            args=[
                "bcftools",
                "view",
                "--no-version",
                "--output-type",
                "z",
                "--output",
                self.file_path,
                "--write-index",
            ],
            stdin=subprocess.PIPE,
            error_message="Creating deidentified records failed",
        )
        self.started = True

    def ingest(self, lines):
        if not self.started:
            self._start()
        print("\n".join(lines), file=self.view_process.stdin)

    def close(self):
        if self.started:
            # self.view_process.stdin.close()
            self.view_process.check()


def anonymise(input_string):
    return ANY_PII_PATTERN.sub(MASK, input_string)


def get_structured_meta_values(value):
    if not (value.startswith("<") and value.endswith(">")):
        raise ParsingError(f"Meta information line is formatted correctly:\n{value}")
    values = {}
    current_key = []
    current_value = []
    extending = current_key
    escaped = False
    in_quotes = False
    for c in value[1:]:
        if escaped:
            extending.append(c)
            escaped = False
        elif c == "\\":
            extending.append(c)
            escaped = True
        elif c == '"':
            in_quotes = not in_quotes
            extending.append(c)
        elif in_quotes:
            extending.append(c)
        elif c in ",>":
            values["".join(current_key)] = "".join(current_value)
            current_key.clear()
            current_value.clear()
            extending = current_key
        elif c == "=":
            extending = current_value
        else:
            extending.append(c)
    return values


def anonymise_header_line(header_line):
    if header_line.startswith("##") and header_line.count("="):
        # Is a meta line
        key, value = header_line[2:].split("=", 1)
        if value.startswith("<"):
            # Structured meta line
            subkey_values = get_structured_meta_values(value)
            if key in META_STRUCTURED_WHITELIST and "Description" in subkey_values:
                subkey_values["Description"] = anonymise(subkey_values["Description"])
            else:
                subkey_values = {
                    anonymise(subkey): anonymise(subvalue)
                    for subkey, subvalue in subkey_values.items()
                }
            new_value = (
                "<"
                + ",".join(
                    f"{subkey}={subvalue}" for subkey, subvalue in subkey_values.items()
                )
                + ">"
            )
            new_line = f"##{key}={new_value}"
        elif key in META_UNSTRUCTURED_WHITELIST:
            new_line = header_line
        else:
            new_line = f"##{anonymise(key)}={anonymise(value)}"
    else:
        # Other comment line or incorrectly formatted, anonymise the whole thing
        new_line = f"#{anonymise(header_line[1:])}"
    return new_line


def anonymise_vcf_record(record, info_whitelist):
    """Anonymise the INFO column of a VCF record"""
    fields = record.split("\t")
    info_field = fields[7]
    info_fields = info_field.split(";")
    new_info_pairs = []
    for field in info_fields:
        if (key_value := field.split("=", 1))[0] not in info_whitelist:
            if len(key_value) == 1:
                # This is a flag field, it should already be in the whitelist
                info_whitelist.add(key_value[0])
            else:
                value = key_value[1]
                new_value = anonymise(value)
                if new_value != value:
                    new_info_pairs.append(f"{key_value[0]}={new_value}")
    if new_info_pairs:
        fields[7] = ";".join(new_info_pairs)
        return "\t".join(fields)
    else:
        return None


def get_output_type(file_path):
    output_type_list = [
        (suffix, output_type)
        for suffix, output_type in GENOMIC_SUFFIX_TYPES.items()
        if file_path.endswith(suffix)
    ]
    assert (
        len(output_type_list) == 1
    ), f"File path {file_path} does not have a valid suffix"
    return output_type_list[0]


def process_header(file_path):
    view_process = CheckedProcess(
        args=["bcftools", "view", "--header-only", "--no-version", file_path],
        stdout=subprocess.PIPE,
        error_message="Reading header failed",
    )
    header_changes = False
    info_whitelist = INFO_RESERVED_KEYS.copy()
    header_lines = []
    for line in view_process.stdout:
        line = line.rstrip("\r\n")
        if line.startswith("##INFO=<"):
            # INFO line, add to whitelist if Type is not "String"
            info_attributes = get_structured_meta_values(line[7:])
            if info_attributes.get("Type", "String") != "String":
                info_whitelist.add(info_attributes.get("ID"))
        new_line = anonymise_header_line(line)
        header_lines.append(new_line)
        if new_line != line:
            header_changes = True
    view_process.check()
    if header_changes:
        print("Header PII detected, creating anonymised header")
        with open(f"{WORKING_DIR}/{HEADER_PATH}", "w") as header_file:
            print("\n".join(header_lines), file=header_file)
    else:
        print("No PII detected in header")
    return info_whitelist, header_lines, header_changes


def process_records(file_path, header_lines, info_whitelist):
    view_process = CheckedProcess(
        args=["bcftools", "view", "--drop-genotypes", "--no-header", file_path],
        stdout=subprocess.PIPE,
        error_message="Reading records failed",
    )
    changed_lines = header_lines.copy()
    # Remove sample columns from header
    changed_lines[-1] = "\t".join(changed_lines[-1].split("\t", 8)[:8])
    num_records_changed = 0
    viewer = Viewer(ANNOTATION_PATH)
    for line in view_process.stdout:
        line = line.rstrip("\r\n")
        new_line = anonymise_vcf_record(line, info_whitelist)
        if new_line is not None:
            changed_lines.append(new_line)
            num_records_changed += 1
            if len(changed_lines) > MAX_LINES_PER_PRINT:
                viewer.ingest(changed_lines)
                changed_lines.clear()
    view_process.check()
    if changed_lines:
        viewer.ingest(changed_lines)
    viewer.close()
    if num_records_changed:
        print(
            f"INFO PII detected, anonymised annotation created for {num_records_changed} record(s)"
        )
    else:
        print("No PII detected in records' INFO columns")
    return num_records_changed > 0


def prepare_for_annotate(file_path):
    """Annotate is very picky, and needs a gzipped indexed file to work"""
    print("Bgzipping and indexing locally for annotation")
    view_process = CheckedProcess(
        args=[
            "bcftools",
            "view",
            "--no-version",
            "--output-type",
            "b",
            "--output",
            BGZIPPED_PATH,
            "--write-index",
            file_path,
        ],
        stdout=subprocess.PIPE,
        error_message="Bgzipping and indexing original file failed",
    )
    view_process.check()


def anonymise_vcf(input_path, output_path):
    output_type = get_output_type(input_path)[1]
    info_whitelist, header_lines, header_changes = process_header(input_path)
    info_changes = process_records(input_path, header_lines, info_whitelist)
    base_reheader_args = [
        "bcftools",
        "reheader",
        "--header",
        HEADER_PATH,
        "--output",
        output_path,
    ]
    base_annotate_args = [
        "bcftools",
        "annotate",
        "--no-version",
        "--annotations",
        ANNOTATION_PATH,
        "--columns",
        "INFO",
        "--pair-logic",
        "exact",
        "--output-type",
        output_type,
        BGZIPPED_PATH,
    ]
    files_to_move = [output_path]
    if output_type in "zb":
        files_to_move.append(f"{output_path}.csi")
    if header_changes:
        if info_changes:
            prepare_for_annotate(input_path)
            reheader_process = CheckedProcess(
                args=base_reheader_args,
                stdin=subprocess.PIPE,
                error_message="Updating header failed",
            )
            annotate_process = CheckedProcess(
                args=base_annotate_args,
                stdout=reheader_process.stdin,
                error_message="Updating INFO column failed",
            )
            reheader_process.check()
            annotate_process.check()
        else:
            reheader_process = CheckedProcess(
                args=base_reheader_args + [input_path],
                error_message="Updating header failed",
            )
            reheader_process.check()
        if output_type in "zb":
            index_process = CheckedProcess(
                args=["bcftools", "index", output_path],
                error_message="Indexing anonymised file failed",
            )
            index_process.check()
    elif info_changes:
        prepare_for_annotate(input_path)
        annotate_process = CheckedProcess(
            args=base_annotate_args
            + ["--output", output_path]
            + (["--write-index"] if output_type in "zb" else []),
            error_message="Updating INFO column failed",
        )
        annotate_process.check()
    else:
        print("No PII detected in VCF file, copying verbatim")
        files_to_move = [input_path]
    return files_to_move


def process_tabular(input_path, output_path, delimiter):
    """Processes CSV/TSV files to deidentify PII, writing results line-by-line."""
    with open(input_path, "r", newline="", encoding="utf-8") as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        with open(output_path, "w", newline="", encoding="utf-8") as outfile:
            writer = csv.writer(outfile, delimiter=delimiter)
            for row in reader:
                deidentified_row = [anonymise(field) for field in row]
                writer.writerow(deidentified_row)


def process_flatfile(input_path, output_path):
    """Processes flat files (.txt, .json) to deidentify PII, writing results line-by-line."""
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            deidentified_line = anonymise(line)
            outfile.write(deidentified_line)


def deidentify_metadata(local_input_path, local_output_path):
    """Main function to process file from S3, deidentify, and upload back to S3."""

    if local_input_path.endswith(".json") or local_input_path.endswith(".txt"):
        process_flatfile(local_input_path, local_output_path)
    elif local_input_path.endswith(".csv"):
        process_tabular(local_input_path, local_output_path, delimiter=",")
    elif local_input_path.endswith(".tsv"):
        process_tabular(local_input_path, local_output_path, delimiter="\t")

    return True


def deidentify(file_name):
    local_input_path = f"{WORKING_DIR}/{file_name}"
    local_output_path = f"{WORKING_DIR}/deidentified_{file_name}"
    if any(file_name.endswith(suffix) for suffix in GENOMIC_SUFFIX_TYPES.keys()):
        try:
            anonymise_vcf(local_input_path, local_output_path)
        except (BcftoolsError, ParsingError) as e:
            print(f"An error occurred while deidentifying {file_name}: {e}")
            print("Exiting")
            return
    elif any(file_name.endswith(suffix) for suffix in METADATA_SUFFIXES):
        deidentify_metadata(local_input_path, local_output_path)
    else:
        raise ValueError(f"File {file_name} does not have a recognised suffix")
    print(f"Successfully deidentify {file_name}")
    print(f"Please Check deidentified_{file_name}")
