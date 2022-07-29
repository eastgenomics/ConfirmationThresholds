from FPTP import *
import unittest
from io import StringIO


class TestRevComp(unittest.TestCase):
    def test_parseArgs(self):
        # not sure how to test this one
        pass


    def test_checkHappyQueryMatch(self):
        self.assertEqual(checkHappyQueryMatch("NA12878.vcf.gz", "NA12878-NA12878-1-TWE-F-EGG4_markdup_recalibrated_Haplotyper.vcf.gz"),True)
        self.assertRaises(AssertionError, checkHappyQueryMatch, "NA73652.vcf.gz" "NA12878-NA12878-1-TWE-F-EGG4_markdup_recalibrated_Haplotyper.vcf.gz")


    def test_checkMultipleQueryMetrics(self):
        self.assertEqual(revComp("accgttaattgccgt"),"acggcaattaacggt")
        self.assertEqual(revComp("ACCGTTAATTGCCGT",True),"ACGGCAATTAACGGT")
        self.assertRaises(SystemExit, revComp, "agtcgahgattc")
        self.assertEqual(revComp("agtcgtagcnnn---taagct"),"agctta---nnngctacgact")


    def test_checkMetrics(self):
        self.assertEqual(revComp("accgttaattgccgt"),"acggcaattaacggt")
        self.assertEqual(revComp("ACCGTTAATTGCCGT",True),"ACGGCAATTAACGGT")
        self.assertRaises(SystemExit, revComp, "agtcgahgattc")
        self.assertEqual(revComp("agtcgtagcnnn---taagct"),"agctta---nnngctacgact")


    def test_parseQuery(self):
        self.assertEqual(revComp("accgttaattgccgt"),"acggcaattaacggt")
        self.assertEqual(revComp("ACCGTTAATTGCCGT",True),"ACGGCAATTAACGGT")
        self.assertRaises(SystemExit, revComp, "agtcgahgattc")
        self.assertEqual(revComp("agtcgtagcnnn---taagct"),"agctta---nnngctacgact")


    def test_parseHappy(self):
        self.assertEqual(revComp("accgttaattgccgt"),"acggcaattaacggt")
        self.assertEqual(revComp("ACCGTTAATTGCCGT",True),"ACGGCAATTAACGGT")
        self.assertRaises(SystemExit, revComp, "agtcgahgattc")
        self.assertEqual(revComp("agtcgtagcnnn---taagct"),"agctta---nnngctacgact")


    def test_getOutputName(self):
        self.assertEqual(revComp("accgttaattgccgt"),"acggcaattaacggt")
        self.assertEqual(revComp("ACCGTTAATTGCCGT",True),"ACGGCAATTAACGGT")
        self.assertRaises(SystemExit, revComp, "agtcgahgattc")
        self.assertEqual(revComp("agtcgtagcnnn---taagct"),"agctta---nnngctacgact")


    def test_makeReport(self):
        pass


    def test_createPlot(self):
        pass


    def test_mergeSamples(self):
        self.assertTrue(checkformat(True,"testFiles/fasta.fasta"))
        self.assertTrue(checkformat(False,"aagggttac"))
        self.assertRaises(AssertionError, checkformat, True,"attaggsc")
        self.assertRaises(AssertionError, checkformat, False,"testFiles/fasta.fasta")
        self.assertFalse(checkformat(True,"testFiles/notfasta.csv"))


    def test_makeArrays(self):
        pass


if __name__ == '__main__':
    unittest.main()