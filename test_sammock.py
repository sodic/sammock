import unittest
from sammock import *



class TestCigarFunctions(unittest.TestCase):
    def test_mark_operation(self):
        self.assertEqual("M", mark_operation("A", "A"))
        self.assertEqual("M", mark_operation("A", "C"))
        self.assertEqual("I", mark_operation("A", "-"))
        self.assertEqual("D", mark_operation("-", "G"))

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

