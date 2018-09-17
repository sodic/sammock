from argparse import ArgumentParser
from itertools import dropwhile, groupby

BLANK_POSITION = "-"
DEFAULT_QUALITY = 60
DEFAULT_REF_NAME = "ref"
VALID_CHARACTERS = set("ACGT").union(BLANK_POSITION)
SAM_ROW_DELIMITER = "\n"
SAM_COL_DELIMITER = "\t"
SAM_VERSION_NUMBER = "1.4"


def is_valid(char):
    return char in VALID_CHARACTERS


def present(base):
    return base != BLANK_POSITION


def missing(base):
    return not present(base)


def remove_all_blanks(sequence):
    return sequence.replace(BLANK_POSITION, "")


def remove_whitespaces(lines):
    return ["".join(line.split()) for line in lines if line.strip()]


def remove_trailing_blanks(sequence):
    return sequence.rstrip(BLANK_POSITION)


def check_legal_characters(lines):
    for row, line in enumerate(lines, 1):
        for col, char in enumerate(line, 1):
            assert is_valid(char), \
                f"Found illegal character '{char}' (line: {row}, position: {col})"


def perform_checks(lines):
    assert len(lines) >= 2, "Provide at least one read and a reference"
    check_legal_characters(lines)


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


def make_sam_entry(index, read, reference):
    shifted_read = "".join(dropwhile(missing, read))
    read_offset = len(read) - len(shifted_read)
    shifted_reference = reference[read_offset:]
    real_read = remove_all_blanks(shifted_read)

    return (index, 0, "ref",
            read_offset + 1, DEFAULT_QUALITY,
            cigar(shifted_read, shifted_reference),
            "*", 0, 0, real_read, "@"*len(real_read))


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


def create_files(reference, reads, ref_file_name, sam_file_name):
    real_reference = remove_all_blanks(reference)
    make_reference_file(real_reference, ref_file_name)

    sam_entries = [make_sam_entry(idx + 1, read, reference)
                   for idx, read in enumerate(reads)]
    sam_entries.sort(key=sam_sorter)

    sam_rows = SAM_COL_DELIMITER.join(str(entry) for entry in sam_entries)
    sam_body = SAM_ROW_DELIMITER.join(sam_rows)

    make_alignments_file(sam_body, len(real_reference), sam_file_name)


def get_args():
    parser = ArgumentParser(description="""Convert a symbolic alignment file 
        into a SAM alignemnt file and a corresponding reference.""")

    parser.add_argument("input_file",
                        help="path to the input file with symbolic alignments.")
    parser.add_argument("-r", "--reference",
                        default="ref.fa",
                        dest="ref_file_name",
                        help="desired name of the generated reference sequence file.")
    parser.add_argument("-a", "--alignments",
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

    lines = remove_whitespaces(lines)
    sequence_strings = [remove_trailing_blanks(l) for l in lines]
    check_legal_characters(sequence_strings)

    reference = sequence_strings.pop()
    reads = sequence_strings

    create_files(reference, reads, ref_file_name, sam_file_name)


if __name__ == "__main__":
    main()
