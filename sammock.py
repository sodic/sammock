import sys
import argparse
from itertools import dropwhile, groupby


def remove_blanks(seq):
    return "".join(c for c in seq if c != "-")


def cigar(read, ref):
    def transform(read_base, ref_base):
        if read_base == '-' and ref_base == '-':
            return ""
        elif read_base == "-" and ref_base != "-":
            return "D"
        elif read_base != "-" and ref_base == "-":
            return "I"
        else:
            return "M"

    grouper = groupby(transform(b1, b2) for b1, b2 in zip(read, ref))
    return "".join(str(len(list(v))) + k for k, v in grouper if k)


def make_sam_entry(index, read, reference):
    shifted_read = "".join(dropwhile(lambda c: c == "-", read))
    read_offset = len(read) - len(shifted_read)
    shifted_reference = reference[read_offset:]
    real_read = remove_blanks(shifted_read)

    return (index, 0, "ref",
            read_offset + 1, 60,
            cigar(shifted_read, shifted_reference),
            "*", 0, 0, real_read, "@"*len(real_read))

def sam_sort_key(entry):
    return entry[2], entry[3]


def format_line(idx, read, reference):
    shifted_read = "".join(dropwhile(lambda c: c == "-", read))
    read_offset = len(read) - len(shifted_read)
    shifted_reference = reference[read_offset:]
    real_read = remove_blanks(shifted_read)

    return "\t".join((f"read{idx}", "0",
                      "ref", str(read_offset + 1), "60",
                      cigar(shifted_read, shifted_reference), "*", "0",
                      "0", real_read, "@"*len(real_read)))


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input_file',
                        help='path to the input file')
    parser.add_argument("-r", "--reference", default="ref.fa", dest="ref_file_name",
                        help="name of the reference sequence file")
    parser.add_argument("-s", "--sam", default="als.sam", dest="sam_file_name",
                        help="name of the alignments (sam) file")

    args = parser.parse_args()
    input_file_name = args.input_file
    ref_file_name = args.ref_file_name
    sam_file_name = args.sam_file_name

    with open(input_file_name) as file:
        content = file.readlines()

    lines = ["".join(line.split()) for line in content if line.strip()]
    assert len(lines) >= 2, "Provide at least one reference and a read"

    valid_chars = set("ACGT-")
    assert all(all(c in valid_chars for c in line)
               for line in lines), "Invalid characters in data"

    reference = lines.pop()
    real_reference = remove_blanks(reference)
    with open(ref_file_name, "w") as ref_file:
        ref_file.write(">ref\n")
        ref_file.write(real_reference)
        ref_file.write("\n")

    SAM_HEADER = f"@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:ref\tLN:{len(real_reference)}\n"
    sam_entries = sorted((make_sam_entry(idx + 1, read, reference) 
            for idx, read in enumerate(lines)), key=sam_sort_key)
    sam_strings = "\n".join("\t".join(map(str,entry)) for entry in sam_entries)
    with open(sam_file_name, "w") as sam_file:
        sam_file.write(SAM_HEADER)
        sam_file.write(sam_strings)


main()
