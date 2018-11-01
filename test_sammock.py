import unittest
from sammock import *
import os


class TestCigarFunctions(unittest.TestCase):
    def test_mark_operation(self):
        self.assertEqual("M", mark_operation("A", "A"))
        self.assertEqual("M", mark_operation("A", "C"))
        self.assertEqual("I", mark_operation("A", "-"))
        self.assertEqual("D", mark_operation("-", "G"))
        self.assertEqual(None, mark_operation("-", "-"))

    def test_cigar_mathces(self):
        self.assertEqual("5M", cigar("ACGCT", "ACGGT"))
        self.assertEqual("3M", cigar("ACT", "ACGGT"))
        self.assertEqual("2M", cigar("ACTTT", "AC"))
        self.assertEqual("1M", cigar("A", "ACGGT"))
        self.assertEqual("1M", cigar("ACCTG", "C"))
        self.assertEqual("1M4D", cigar("A----", "ACGGT"))
        self.assertEqual("1M1I", cigar("AD", "A----"))

    def test_cigar_insertions(self):
        self.assertEqual("1M2D2M", cigar("A--CT", "ACGGT"))
        self.assertEqual("1M3D1M", cigar("A---T", "ACGGT"))
        self.assertEqual("2I2M", cigar("CGTT", "--TTT"))

    def test_cigar_deletions(self):
        self.assertEqual("1M2I2M", cigar("ADGCG", "A--CG"))
        self.assertEqual("1M3I1M", cigar("ADGCG", "A---G"))
        self.assertEqual("2D2M", cigar("--TC", "ACGGT"))

    def test_cigar_edge(self):
        self.assertEqual("", cigar("", "CGC--TA---AGATTGT"))
        self.assertEqual("", cigar("GCGGD", ""))

    def test_cigar_complex(self):
        self.assertEqual("2M2I1M2D3M", cigar("GCTGA--GTC-", "AC--AGCAGG--GC"))
        self.assertEqual("4D1M2I6M1D1M2D1M", cigar("------CCATA---AGTA-C--T", "--CACGC--TA---AGTATTATTGT"))
        reference = "A-GA-GAC-AC"
        for cg, read in [("1D1I2M1I3M", "-CGTACCT"),
                         ("2M1D3M", "A-G--GTA"),
                         ("1D2M1I3M", "--TGACGT"),
                         ("1M1I2M1I2D1M1I", "CGTAC--AC"),
                         ("1M1I1D1M1I3M1I1M", "AC-TCGATGA"),
                         ("1M1I1M1D2M1D2M", "ACG--GA--AC"),
                         ("1D1I4M1D2M", "-CGA-GA--AC")]:
            self.assertEqual(cg, cigar(read, reference))


class TestStringManipulations(unittest.TestCase):
    def test_remove_whitespaces(self):
        self.assertEqual(remove_whitespaces("A C G T"), "ACGT")
        self.assertEqual(remove_whitespaces("  A\t C    G\n T"), "ACGT")

    def test_parse_base(self):
        self.assertEqual(("A", 34), parse_base("A:34"))
        self.assertEqual(("C", 52), parse_base("C:52"))
        self.assertEqual(("G", None), parse_base("G"))
        self.assertEqual((BLANK_POSITION, None), parse_base(BLANK_POSITION))
        self.assertEqual((BLANK_POSITION, None), parse_base(BLANK_POSITION))

    def test_legal_reads_legals(self):
        self.assertEqual(True,
                         is_legal_read( "A:42 A:42 A:42 A:42 - G:43 C:45 C:42 T:44 T:44 A:43 - - A:42 A:43"))
        self.assertEqual(True, is_legal_read("A C - - G C G - T"))
        self.assertEqual(True, is_legal_read("A:30 C:40 - - G:52 C:21 G:34 - T:36"))
        self.assertEqual(True, is_legal_read("- A:30 C:40 - - G:52 C:21 G:34 - T:36 -"))

    def test_legal_reads_illegal(self):
        self.assertEqual(False, is_legal_read("A C:40 - - G C G - T"))
        self.assertEqual(False, is_legal_read("A C - - G U G - T"))
        self.assertEqual(False, is_legal_read("A C -:32 - G C G - T"))
        self.assertEqual(False, is_legal_read("A:23 C:40 - - G:20 C:21 G - T:43"))

    def test_quality_string(self):
        self.assertEqual("8LUB", quality_string([23, None, 43, 52, 33]))
        self.assertEqual("8NLUB", quality_string([23, 45, 43, 52, 33]))

    def test_parse_read(self):
        result = [("C", None), ("A", 34), ("C", 44)] + \
                 [(BLANK_POSITION, None)] * 3 + \
                 [("G", 65), ("G", None), ("A", 42)] + \
                 [(BLANK_POSITION, None), ("A", None)]
        self.assertEqual(result, parse_read("C A:34 C:44 - - -  G:65 G A:42 - A"))


class TestSammock(unittest.TestCase):
    def test_samples(self):
        current = os.path.dirname(__file__)
        tests = os.path.join(current, "test")

        for sample_dir in os.listdir(tests):
            abosulute = os.path.join(tests, sample_dir)

            fasta_expected = read_file(os.path.join(abosulute, "ref.fa"))
            sam_expected = read_file(os.path.join(abosulute, "alignments.sam"))

            content = read_file(os.path.join(abosulute, "input.txt"))
            fasta_acutal, sam_actual = sammock(content)
            write_file(sam_actual, "test.sam")

            self.assertEqual(fasta_expected, fasta_acutal)
            self.assertEqual(sam_expected, sam_actual)
