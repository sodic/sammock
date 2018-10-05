#!/usr/bin/python3.6
from argparse import ArgumentParser
from itertools import dropwhile, groupby
from re import compile

BLANK_POSITION = "-"
DEFAULT_QUALITY = 60
DEFAULT_REF_NAME = "ref"
VALID_CHARACTERS = "ACGT"
BLANK_POSITION = "-"
SAM_ROW_DELIMITER = "\n"
SAM_COL_DELIMITER = "\t"
SAM_VERSION_NUMBER = "1.4"
LEGAL_READ_PATTERN = compile(r"(([" + VALID_CHARACTERS + "](:[0-9]*){0,1}|\-)([\s]+|$))*")


def present(base):
    return base != BLANK_POSITION


def missing(base):
    return not present(base)


def remove_all_blanks(sequence):
    return sequence.replace(BLANK_POSITION, "")


def remove_whitespaces(line):
    return "".join(line.strip().split()) 


def remove_trailing_blanks(sequence):
    return sequence.rstrip(BLANK_POSITION)


def is_legal_read(read_str):
    return LEGAL_READ_PATTERN.match(read_str)


def check_legal_reads(read_strs):
    for index, read in enumerate(read_strs, 1):
        if not is_legal_read(read):
            raise ValueError(f"Illegal read on line number {index}.")


def read_string(read_data):
    return (base for base, _ in read_data)


def read_quals(read_data):
    return (quality for _, quality in read_data)


def mark_operation(read_base, ref_base):
    if missing(read_base) and missing(ref_base):
        return None
    elif missing(read_base) and present(ref_base):
        return "D"
    elif present(read_base) and missing(ref_base):
        return "I"
    else:
        return "M"


def cigar_entry(operation, iterator):
    return str(sum(1 for _ in iterator)) + operation


def cigar(read, ref):
    operations = (mark_operation(*bases) for bases in zip(read, ref))
    grouper = groupby(op for op in operations if op)
    return "".join(cigar_entry(op, values) for op, values in grouper if op)


def make_sam_entry(index, read_data, reference_string):
    read_string, read_qualities = zip(*read_data)
    read_qualities = [q for q in read_qualities if q]

    shifted_read_string = "".join(dropwhile(missing, read_string))
    read_offset = len(read_string) - len(shifted_read_string)
    real_read_string = remove_all_blanks(shifted_read_string)

    shifted_reference_string = reference_string[read_offset:]

    return (index, 0, "ref",
            read_offset + 1, DEFAULT_QUALITY,
            cigar(shifted_read_string, shifted_reference_string),
            "*", 0, 0, real_read_string, "".join(chr(q + 33) for q in read_qualities))


def make_reference_file(ref_sequence, ref_file_name):
    with open(ref_file_name, "w") as ref_file:
        ref_file.write(f">{DEFAULT_REF_NAME}\n")
        ref_file.write(ref_sequence)
        ref_file.write("\n")


def make_sam_header(ref_length):
    SAM_HEADER_R1 = SAM_COL_DELIMITER.join([
        "@HD",
        f"VN:{SAM_VERSION_NUMBER}",
        "SO:coordinate"
    ])
    SAM_HEADER_R2 = SAM_COL_DELIMITER.join([
        "@SQ",
        f"SN:{DEFAULT_REF_NAME}",
        f"LN:{ref_length}"
    ])
    return SAM_ROW_DELIMITER.join([SAM_HEADER_R1, SAM_HEADER_R2])


def make_alignments_file(sam_body, ref_length, sam_file_name):
    sam_header = make_sam_header(ref_length)

    with open(sam_file_name, "w") as sam_file:
        sam_file.write(sam_header)
        sam_file.write(SAM_ROW_DELIMITER)
        sam_file.write(sam_body)


def sam_sorter(entry):
    return entry[2], entry[3]


def make_base_info(string):
    info = string.split(":")

    base = info[0]
    quality = info[1] if len(info) == 2 else DEFAULT_QUALITY

    return base, int(quality)
    

def create_files(reference, reads, ref_file_name, sam_file_name):
    real_reference = remove_all_blanks(reference)
    make_reference_file(real_reference, ref_file_name)
   
    reads_data = [[make_base_info(base_str) for base_str in read.split()]
            for read in reads]

    sam_entries = [make_sam_entry(idx + 1, read_data, reference)
                   for idx, read_data in enumerate(reads_data)]
    sam_entries.sort(key=sam_sorter)

    sam_body = SAM_ROW_DELIMITER.join(SAM_COL_DELIMITER.join(
        map(str, entry)) for entry in sam_entries)

    make_alignments_file(sam_body, len(real_reference), sam_file_name)


def get_args():
    parser = ArgumentParser(description="""Convert a symbolic alignment file
        into a SAM alignemnt file and a corresponding reference.""")

    parser.add_argument(
        "input_file",
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


def main():
    args = get_args()
    input_file_name = args.input_file
    ref_file_name = args.ref_file_name
    sam_file_name = args.sam_file_name

    with open(input_file_name) as file:
        lines = file.readlines()
    
    assert len(lines) >= 2, "Provide at least one read and a reference"

    sequence_strings = [line.strip() for line in lines]

    reference_string = remove_trailing_blanks(
            remove_whitespaces(sequence_strings.pop()))
    
    reads = sequence_strings
    check_legal_reads(reads)

    create_files(reference_string, reads, ref_file_name, sam_file_name)


if __name__ == "__main__":
    main()
