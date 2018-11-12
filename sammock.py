#!/usr/bin/python3.6
from argparse import ArgumentParser
from itertools import dropwhile, groupby
from re import compile
from typing import Pattern, List, Iterable, Tuple, Optional

BaseInfo = Tuple[str, Optional[int]]
SamEntry = Tuple[int, int, str, int, int, str, str, int, int, str, str]

BLANK_POSITION: str = "-"
DEFAULT_QUALITY: int = 60
DEFAULT_REF_NAME: str = "ref"

SAM_ROW_DELIMITER: str = "\n"
SAM_COL_DELIMITER: str = "\t"
SAM_VERSION_NUMBER: str = "1.4"

VALID_CHARACTERS: str = "ACGT"
SPACES_OR_LINE_END = r"[\s]+|$"
READ_WITH_QUALITIES = fr"(([{VALID_CHARACTERS}]:[0-9]+|-)({SPACES_OR_LINE_END}))+\Z"
READ_WITHOUT_QUALITIES = fr"(([{VALID_CHARACTERS}]|-)({SPACES_OR_LINE_END}))+\Z"
LEGAL_READ = f"({READ_WITH_QUALITIES})|({READ_WITHOUT_QUALITIES})"
LEGAL_READ_PATTERN: Pattern[str] = compile(LEGAL_READ)


def present(base: str) -> bool:
    return base != BLANK_POSITION


def missing(base: str) -> bool:
    return not present(base)


def remove_all_blanks(sequence: str) -> str:
    return sequence.replace(BLANK_POSITION, "")


def remove_whitespaces(line: str) -> str:
    return "".join(line.strip().split())


def remove_trailing_blanks(sequence: str) -> str:
    return sequence.rstrip(BLANK_POSITION)


def is_legal_read(read: str) -> bool:
    return bool(LEGAL_READ_PATTERN.match(read))


def check_legal_reads(read_strs: Iterable[str]) -> None:
    for index, read in enumerate(read_strs, 1):
        if not is_legal_read(read):
            raise ValueError(f"Illegal read on line number {index}.")


def mark_operation(read_base: str, ref_base: str) -> Optional[str]:
    if missing(read_base) and missing(ref_base):
        return None
    elif missing(read_base) and present(ref_base):
        return "D"
    elif present(read_base) and missing(ref_base):
        return "I"
    else:
        return "M"


def cigar_entry(operation: str, iterator: Iterable) -> str:
    return str(sum(1 for _ in iterator)) + operation


def cigar(read: str, ref: str) -> str:
    operations = (mark_operation(*bases) for bases in zip(read, ref))
    grouper = groupby(op for op in operations if op)
    return "".join(cigar_entry(op, values) for op, values in grouper if op)


def quality_string(qualities: List[int]) -> str:
    result = "".join(chr(q + 33) for q in qualities if q is not None)
    return result if result else "*"


def value_and_qual_strings(read: List[BaseInfo]) -> Tuple[str, str]:
    read_str, read_qualities = zip(*read)
    return "".join(read_str).rstrip(BLANK_POSITION), \
           quality_string(read_qualities)


def adjust_for_reference_blanks(read_offset: int, reference: str) -> int:
    return read_offset - reference[:read_offset].count(BLANK_POSITION)


def make_sam_entry(index: int,
                   read: List[BaseInfo],
                   reference: str) -> SamEntry:
    read_str, qual_str = value_and_qual_strings(read)

    shifted_read = "".join(dropwhile(missing, read_str))
    real_read = remove_all_blanks(shifted_read)
    read_offset = len(read_str) - len(shifted_read)

    shifted_reference = reference[read_offset:]

    real_offset = adjust_for_reference_blanks(read_offset, reference)

    return (index, 0, "ref", real_offset + 1, DEFAULT_QUALITY,
            cigar(shifted_read, shifted_reference),
            "*", 0, 0, real_read, qual_str)


def write_file(content: str, file_name: str) -> None:
    with open(file_name, "w") as file:
        file.write(content)


def make_sam_header(ref_length: int) -> str:
    sam_header_r1 = SAM_COL_DELIMITER.join([
        "@HD",
        f"VN:{SAM_VERSION_NUMBER}",
        "SO:coordinate"
    ])
    sam_header_r2 = SAM_COL_DELIMITER.join([
        "@SQ",
        f"SN:{DEFAULT_REF_NAME}",
        f"LN:{ref_length}"
    ])
    return SAM_ROW_DELIMITER.join([sam_header_r1, sam_header_r2])


def sam_sorter(entry: SamEntry) -> Tuple[str, int]:
    return entry[2], entry[3]


def parse_base(string: str) -> BaseInfo:
    info = string.split(":")

    base = info[0]
    quality = int(info[1]) if len(info) == 2 else None

    return base, quality


def parse_read(read: str) -> List[BaseInfo]:
    return [parse_base(base_str) for base_str in read.split()]


def make_sam_body(entries: List[SamEntry]) -> str:
    return SAM_ROW_DELIMITER.join(SAM_COL_DELIMITER.join(
        map(str, entry)) for entry in entries)


def make_fasta_string(reference: str) -> str:
    return f">{DEFAULT_REF_NAME}\n{reference}\n"


def make_sam_string(reference: str, reads: List[str]) -> str:
    real_reference = remove_all_blanks(reference)
    reads = [parse_read(read) for read in reads]

    sam_entries = [make_sam_entry(idx + 1, read, reference)
                   for idx, read in enumerate(reads)]
    sam_entries.sort(key=sam_sorter)

    sam_body = make_sam_body(sam_entries)
    sam_header = make_sam_header(len(real_reference))

    return SAM_ROW_DELIMITER.join([sam_header, sam_body])


def sammock(input_str: str) -> Tuple[str, str]:
    lines = input_str.split("\n")

    # we are ignoring empty lines
    lines = (line.strip() for line in lines)
    sequences = [line for line in lines if len(line) > 0]

    assert len(sequences) >= 2, "Provide at least one read and a reference"

    reference = remove_trailing_blanks(remove_whitespaces(sequences.pop()))
    real_reference = remove_all_blanks(reference)

    reads = sequences
    check_legal_reads(reads)

    return make_fasta_string(real_reference), make_sam_string(reference, reads)


def read_file(input_file_name):
    with open(input_file_name) as file:
        return file.read()


def get_args():
    parser = ArgumentParser(description="""Convert a symbolic alignment file
        into a SAM alignemnt file and a corresponding reference.""")

    parser.add_argument(
        "input_file_name",
        help="path to the input file with symbolic alignments.")
    parser.add_argument(
        "-r",
        "--reference",
        default="ref.fa",
        dest="ref_file_name",
        help="desired name of the generated reference sequence file.")
    parser.add_argument(
        "-a",
        "--alignments",
        default="alignments.sam",
        dest="sam_file_name",
        help="desired name of the generated alignments (sam) file.")

    return parser.parse_args()


def main(input_file_name, ref_file_name, sam_file_name):
    input_file_content = read_file(input_file_name)

    ref_file_content, sam_file_content = sammock(input_file_content)

    write_file(sam_file_content, sam_file_name)
    write_file(ref_file_content, ref_file_name)


if __name__ == "__main__":
    args = get_args()
    main(args.input_file_name, args.ref_file_name, args.sam_file_name)
