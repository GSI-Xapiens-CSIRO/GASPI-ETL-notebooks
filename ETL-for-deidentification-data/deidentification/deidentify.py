"""
Notebook version of the deidentification pipeline.

Logic ported 1:1 from the application's `deidentification.py`, with the
following intentionally stripped out (not relevant for notebook experiments):
  - S3 download/upload
  - DynamoDB status / error logging
  - CloudWatch marker prints
  - argparse CLI entrypoint
  - BAM/SAM support (notebook input is typically VCF/BCF or small metadata)
  - file_validation (python-magic + htsfile)

Dependencies:
  - bcftools / htslib binary on PATH (for VCF/BCF processing)
  - Python package: ijson  (pip install ijson)

Usage (from the notebook):
    from deidentification.deidentify import deidentify
    deidentify("PII.json")

Output is written next to the input as `deidentified_<filename>`.
"""

import csv
import io
import json
import os
import re
import subprocess

import ijson


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

    # Phone (HP Indonesia)
    r"\b(?:\+?62[\s.-]*|0)8[1-9](?:[\s.-]?\d){7,10}\b",

    # Landline with area code
    r"\b0(?:\((?:2[1-9]|[3-9]\d)\)|(?:2[1-9]|[3-9]\d))(?:[\s.-]?\d){6,8}\b",

    # NIK (without anchor, can be used for SIM also)
    r"\b(1[1-9]|21|[37][1-6]|5[1-3]|6[1-5]|[89][12])\d{2}\d{2}([04][1-9]|[1256][0-9]|[37][01])(0[1-9]|1[0-2])\d{2}\d{4}\b",

    r"\b[A-Z]{1,2} \d{1,4}( [A-Z]{1,3})?\b",  # License plate with enforced spaces
    r"\b(?:[A-Z]\d{6,7}[A-Z]?|[A-Z]{2}\d{6,7})\b",  # Passport
    r"\b\d{13}\b",  # BPJS

    r"\b\d{2}\.?\d{3}\.?\d{3}\.?\d{1}-?\d{3}\.?\d{3}\b",  # NPWP 15 digit
    r"\b0\d{15}\b",  # NPWP 16 digit for foreigners

    # Lat/Long pair
    r"(?<!\w)[+-]?(?:90(?:\.0+)?|[0-8]?\d(?:\.\d+)?)(?:\s*,\s*|\s+)[+-]?(?:180(?:\.0+)?|1[0-7]\d(?:\.\d+)?|[0-9]?\d(?:\.\d+)?)(?!\w)",

    # Single lat/long value
    r"(?<!\w)[+-]?(?:90(?:\.0+)?|[0-8]?\d)\.\d{4,}(?!\w)|(?<!\w)[+-]?(?:180(?:\.0+)?|1[0-7]\d|[0-9]?\d)\.\d{4,}(?!\w)",
]
CASE_INSENSITIVE_PII_PATTERNS = [
    r"\b(?:(?:Jl\.|Jalan|Desa|Kelurahan|Kecamatan|Kab.|Kab|Kabupaten|Kec|Kec.|Kecamatan|Prov.|Provinsi|Prov|Kode\s?Pos)(?:\s?(?:\d{5}|RT\s?\d{1,2}/RW\s?\d{1,2}|[A-Z^RT]+[a-z]*(?:\.\s?\d+)?),?)+,?\s?)+\b",  # Address
    r"\b(?:Dr\.|Prof\.|Ir\.|Haji|Hajjah|Putra|Putri|Sri|Adi|Raden|Ny|H\.|Hj\.|Kiai|Kyai|K\.H\.|KH\.|Gus|Ning|Ir\.|Drs\.|Dra\.|Sultan|Pangeran|R\.M\.|R\.A\.|R\.Ay\.|Rr\.|Ust\.|Ustaz|Ustadz|Bapak|Bpk\.|Iibu|Saudara|Saudari|Sdr\.|Tuan|Tn\.|Nona|Nn\.|Dokter|Drg\.|Prof\.|Teungku|Teuku|Tgk\.|Datuk|Datuak|Tengku|Kemas|Nyimas|Kiagus|Nyanyu|Tubagus|Ratu|Rd\.|Rara|Roro|Anak Agung|Gusti|Dewa|Desak|Lalu|Umbu|ADaeng|Puang|Kapitan|Adi|Daeng|Karaeng|Arung|Opu|Petta|Latu|Upu)(?:\s[A-Z][a-z]+){1,2}\b",  # Name
]
ANY_PII_PATTERN = re.compile(
    "|".join(
        f"(?:{pattern})"
        for pattern in PII_PATTERNS
        + [f"(?i:{pattern})" for pattern in CASE_INSENSITIVE_PII_PATTERNS]
    )
)
INDIVIDUAL_MARKER_FIELDS = {
    "age",
    "umur",
    "ethnicity",
    "etnis",
    "karyotypicsex",
    "sex",
    "jenis kelamin",
    "jenis_kelamin",
}  # What we use to know we can use NAME_PATTERN to remove any name fields
NAME_PATTERN = re.compile(
    r"(?i)(?:^|[_\s])(?:name|nama|marga|initial|inisial|"
    r"nama_lengkap|fullname|full_name|nama_depan|nama_belakang|first_name|last_name|nama_ibu|nama ibu|nama_ayah|nama ayah|nama_pasangan|nama pasangan|nama_wali|nama wali|wali|kontak_darurat|kontak darurat|emergency_contact|gelar|title)(?:$|[_\s])"
)  # Very broad - need to know we're dealing with an individual to use this
METADATA_KEY_PII_PATTERNS = [
    r"(?i)\b(?:(?:full|first|last|middle|given|family|sur)[_ -]?name|nama(?:[_ -](?:lengkap|depan|belakang|tengah))?|nama|surname)\b",
    r"(?i)\b(?:(?:plate|license|vehicle|registration|number)_(?:plate|number|nopol|polisi|registrasi)|(?:nomor|plat)_(?:plat|nomor|polisi|registrasi)|nopol(?:_id)?|vehicle_nopol|registration_nopol|plat_number|plateno)\b",
    r"(?i)\b(?:dob|date[_ -]*of[_ -]*birth|birth[_ -]*date|birthdate|tanggal[_ -]*(?:lahir|lhr)|tgl[_ -]*(?:lahir|lhr))\b",
    r"(?i)\b(?:birth[_ -]*place|place[_ -]*of[_ -]*birth|tempat[_ -]*(?:lahir|lhr)|tmp[_ -]*(?:lahir|lhr))\b",
    r"(?i)\b(?:lokasi)\b",
    r"(?i)\b(?:alamat(?:[_ -]*(?:lengkap|rumah|domisili|ktp|tempat[_ -]*tinggal|tinggal))?|tempat[_ -]*tinggal|tinggal|domisili|rumah|address|full[_ -]*address|home[_ -]*address|domicile|home)\b",
    r"(?i)\b(?:gps|koordinat|coordinate|coordinates|latitude|lat|longitude|lon|lng)\b",
    r"(?i)\b(?:agama|religion|kepercayaan)\b",
    r"(?i)\b(?:(?:no|nomor|nomer)[_ -]*(?:telepon|telfon|telphone|telephone|phone|hp|handphone|ponsel|mobile|whatsapp|wa)|(?:telepon|telfon|telphone|telephone|phone|hp|handphone|ponsel|mobile|whatsapp|wa))\b",
    r"(?i)\b(?:(?:no|nomor|nomer)[_ -]*sim|sim)\b",
    r"(?i)\b(?:mrn|medical[_ -]*record[_ -]*number|nomor[_ -]*rekap[_ -]*medis|no[_ -]*rm|nomor[_ -]*rm|no[_ -]*mr|nomor[_ -]*mr|mr[_ -]*number)\b",
]

CONTROLLED_TERM_PREFIXES = (
    "ICD9CM", "ICD10", "LOINC", "SNOMED", "UCUM", "KFA", "NCIT", "KEMKES"
)

CONTROLLED_TERM_PATTERN = re.compile(
    rf"^(?:{'|'.join(CONTROLLED_TERM_PREFIXES)}):\S+$"
)

GENOMIC_SUFFIX_TYPES = {
    ".bcf": "u",
    ".bcf.gz": "b",
    ".bcf.bgz": "b",
    ".vcf": "v",
    ".vcf.gz": "z",
    ".vcf.bgz": "z",
}

METADATA_SUFFIXES = [
    ".json",
    ".csv",
    ".tsv",
    ".txt",
]

WORKING_DIR = os.getcwd()


class ProcessError(Exception):
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
            raise ProcessError(
                self.error_message, stdout, stderr, returncode, self.process.args
            )


class Viewer:
    def __init__(self, process_args, error_message):
        self.started = False
        self.process_args = process_args
        self.error_message = error_message
        self.view_process = None
        self.lines = []

    def _print(self, lines):
        if not self.started:
            self._start()
        try:
            print("\n".join(lines), file=self.view_process.stdin)
        except BrokenPipeError:
            self.view_process.check()
            # If that's not the cause, raise for further inspection
            raise

    def _start(self):
        self.view_process = CheckedProcess(
            args=self.process_args,
            stdin=subprocess.PIPE,
            error_message=self.error_message,
        )
        self.started = True

    def ingest(self, new_lines):
        self.lines.extend(new_lines)
        if len(self.lines) > MAX_LINES_PER_PRINT:
            self._print(self.lines)
            self.lines.clear()

    def close(self):
        if self.lines:
            self._print(self.lines)
        if self.started:
            self.view_process.check()


def anonymise(input_string):
    if CONTROLLED_TERM_PATTERN.match(input_string):
        return input_string
    return ANY_PII_PATTERN.sub(MASK, input_string)


def remove_nested_angle_brackets(header_line):
    match = re.search(r"^(##\w+)=<(.+)>$", header_line)
    if not match:
        return header_line
    prefix, inner = match.groups()
    flattened_inner = inner.replace("<", "").replace(">", "")
    return f"{prefix}=<{flattened_inner}>"


def get_structured_meta_values(value):
    if not (value.startswith("<") and value.endswith(">")):
        raise ParsingError(f"Meta information line is formatted incorrectly:\n{value}")
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
        full_length = len(line)
        line = line.rstrip("\r\n")
        if len(line) == full_length:
            # No line ending, has view_process crashed?
            view_process.check()
        line = line.rstrip("\r\n")
        if line.startswith("##INFO=<"):
            # INFO line, add to whitelist if Type is not "String"
            info_attributes = get_structured_meta_values(line[7:])
            if info_attributes.get("Type", "String") != "String":
                info_whitelist.add(info_attributes.get("ID"))
        new_line = anonymise_header_line(remove_nested_angle_brackets(line))
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
    header_lines = header_lines.copy()
    # Remove sample columns from header
    header_lines[-1] = "\t".join(header_lines[-1].split("\t", 8)[:8])
    num_records_changed = 0
    viewer = Viewer(
        [
            "bcftools",
            "view",
            "--no-version",
            "--output-type",
            "z",
            "--output",
            ANNOTATION_PATH,
            "--write-index",
        ],
        "Creating deidentified records failed",
    )
    for line in view_process.stdout:
        full_length = len(line)
        line = line.rstrip("\r\n")
        if len(line) == full_length:
            # No line ending, has view_process crashed?
            view_process.check()
        line = line.rstrip("\r\n")
        new_line = anonymise_vcf_record(line, info_whitelist)
        if new_line is not None:
            if num_records_changed == 0:
                viewer.ingest(header_lines + [new_line])
            else:
                viewer.ingest([new_line])
            num_records_changed += 1
    view_process.check()
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
                args=["bcftools", "index", "--force", output_path],
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
        if output_type in "zb":
            index_process = CheckedProcess(
                args=["bcftools", "index", "--force", input_path],
                error_message="Indexing original file failed",
            )
            index_process.check()
            files_to_move.append(f"{input_path}.csi")
    return files_to_move


def process_tabular(input_path, output_path, delimiter):
    """Processes CSV/TSV files to deidentify PII and drop sensitive columns."""
    with open(input_path, "r", newline="", encoding="utf-8") as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        header = next(reader)
        is_individual = any(
            col_name.casefold() in INDIVIDUAL_MARKER_FIELDS for col_name in header
        )
        columns_to_keep = [
            idx
            for idx, col_name in enumerate(header)
            if not any(
                re.match(pattern, col_name) for pattern in METADATA_KEY_PII_PATTERNS
            )
            and not (is_individual and NAME_PATTERN.search(col_name))
        ]
        with open(output_path, "w", newline="", encoding="utf-8") as outfile:
            writer = csv.writer(outfile, delimiter=delimiter)

            filtered_header = [header[idx] for idx in columns_to_keep]
            writer.writerow(filtered_header)
            for row in reader:
                filtered_row = [anonymise(row[idx]) for idx in columns_to_keep]
                writer.writerow(filtered_row)


def outfile_after_element(stack, outfile):
    if stack:
        if "stored_outfile" in stack[-1]:
            stack[-1].setdefault("name_strings", []).append(
                outfile.getvalue().strip(",")
            )
            outfile.close()
            outfile = stack[-1].pop("stored_outfile")
        if stack[-1].get("skip_first"):
            del stack[-1]["skip_first"]
        else:
            stack[-1]["first"] = False
    return outfile


def process_json(input_path, output_path):
    """Process JSON files to deidentify PII, writing results line-by-line and omitting sensitive keys"""
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        parser = ijson.parse(infile)
        # The stack holds a dictionary for each container with keys:
        #  'type': "object" or "array"
        #  'first': boolean flag, True if no item has been written yet.
        #  'pending_key': for objects, True if a key was written but its value not has not yet been written.
        stack = []
        keybuffer = None  # When set, skip all subelements

        for prefix, event, value in parser:
            if keybuffer and not prefix.startswith(keybuffer):
                keybuffer = None
            if keybuffer:
                continue

            if event == "start_map":
                if stack:
                    if stack[-1]["type"] == "object" and stack[-1].get("pending_key"):
                        outfile.write(":")
                        stack[-1]["pending_key"] = False
                    elif stack[-1]["type"] == "array" and not stack[-1]["first"]:
                        outfile.write(",")
                outfile.write("{")
                stack.append({"type": "object", "first": True, "pending_key": False})

            elif event == "end_map":
                if "name_strings" in stack[-1] and not stack[-1].get("is_individual"):
                    if not stack[-1]["first"]:
                        outfile.write(",")
                    outfile.write(",".join(stack[-1].pop("name_strings")))
                outfile.write("}")
                stack.pop()
                outfile = outfile_after_element(stack, outfile)

            elif event == "start_array":
                if stack:
                    if stack[-1]["type"] == "object" and stack[-1].get("pending_key"):
                        outfile.write(":")
                        stack[-1]["pending_key"] = False
                    elif stack[-1]["type"] == "array" and not stack[-1]["first"]:
                        outfile.write(",")
                outfile.write("[")
                stack.append({"type": "array", "first": True})

            elif event == "end_array":
                outfile.write("]")
                stack.pop()
                outfile = outfile_after_element(stack, outfile)

            elif event == "map_key":
                if value.casefold() in INDIVIDUAL_MARKER_FIELDS:
                    stack[-1]["is_individual"] = True
                # If the key matches a PII pattern, set the keybuffer to skip its subtree.
                if any(
                    re.match(pattern, value) for pattern in METADATA_KEY_PII_PATTERNS
                ):
                    keybuffer = f"{prefix}.{value}"
                    continue
                if NAME_PATTERN.search(value):
                    stack[-1]["stored_outfile"] = outfile
                    outfile = io.StringIO()
                    stack[-1]["skip_first"] = True
                if stack and stack[-1]["type"] == "object":
                    if not stack[-1]["first"]:
                        outfile.write(",")
                    outfile.write(json.dumps(value))
                    stack[-1]["pending_key"] = True

            elif event in ("string", "number", "boolean", "null"):
                if stack:
                    if stack[-1]["type"] == "object" and stack[-1].get("pending_key"):
                        outfile.write(":")
                        stack[-1]["pending_key"] = False
                    elif stack[-1]["type"] == "array":
                        if not stack[-1]["first"]:
                            outfile.write(",")
                if event == "string":
                    outfile.write(json.dumps(anonymise(value)))
                else:
                    outfile.write(json.dumps(value))
                outfile = outfile_after_element(stack, outfile)

        outfile.write("\n")


def process_flatfile(input_path, output_path):
    """Processes TXT files to deidentify PII, writing results line-by-line."""
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            deidentified_line = anonymise(line)
            outfile.write(deidentified_line)


def deidentify_metadata(local_input_path, local_output_path):
    """Dispatch metadata files to the right processor based on suffix."""

    if local_input_path.endswith(".json"):
        process_json(local_input_path, local_output_path)
    elif local_input_path.endswith(".txt"):
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
        except (ProcessError, ParsingError) as e:
            print(f"An error occurred while deidentifying {file_name}: {e}")
            print("Exiting")
            return
    elif any(file_name.endswith(suffix) for suffix in METADATA_SUFFIXES):
        deidentify_metadata(local_input_path, local_output_path)
    else:
        raise ValueError(f"File {file_name} does not have a recognised suffix")
    print(f"Successfully deidentify {file_name}")
    print(f"Please Check deidentified_{file_name}")
