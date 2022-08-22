import FPTP
import unittest
import json
import numpy as np
import os

# these test objects are correct as long as the test input files bundled with
# the module are unchanged
with open('test_data.json', 'r') as fh:
    test_data = json.load(fh)

test_happy_dict = test_data['test_happy_dict']
test_query_dict = test_data['test_query_dict']
test_merged_dict = test_data['test_merged_dict']
test_lists = test_data['test_lists']


class TestModule(unittest.TestCase):
    def test_calculate_centiles(self):
        '''
        Checks that the calculate_centiles() function returns the correct
        output.

        TODO - change expected output to numpy array
        '''
        a = FPTP.calculate_centiles([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        b = np.array(
            [9.09, 18.18, 27.27, 36.36, 45.45, 54.55, 63.64, 72.73, 81.82,
             90.91, 100.]
             )
        self.assertEqual(
            np.array_equal(a, b), True
            )

    def test_check_happy_query_match(self):
        '''
        Checks that the check_happy_query_match() function returns True if the
        sample names match, and raises an AssertionError if not.
        '''
        self.assertEqual(
            FPTP.check_happy_query_match("NA12878", "NA12878"), True
            )
        self.assertRaises(
            AssertionError, FPTP.check_happy_query_match, "NA12878", "sample1"
            )

    def test_check_metrics(self):
        '''
        Checks that the check_metrics() function returns the correct list of
        metrics, either what was requested or this minus metrics not present
        in the query vcf.
        '''
        a = sorted(
            FPTP.check_metrics(
                "NA12878-query-test.vcf.gz", 'AD,DP,GQ,AC,AF,AN,BaseQRankSum'
                )
            )
        b = sorted(
            FPTP.check_metrics(
                "NA12878-query-test.vcf.gz",
                'AD,DP,GQ,AC,AF,AN,BaseQRankSum,FictionalMetric,BDA'
                )
            )
        c = sorted(
            [
                'info_AN', 'info_DP', 'info_BaseQRankSum',
                'format_GQ', 'format_DP'
                ]
            )
        self.assertEqual(a, c)
        self.assertEqual(b, c)

    def test_check_multiple_query_metrics(self):
        # not implemented yet (multiple samples)
        pass

    def test_create_plot(self):
        '''
        Check that create_plot() returns a plotly figure object matching
        the example.
        '''
        # discussed with Matt - can't test something so big as no orthogonal
        # method for making the example
        pass

    def test_decide_bins(self):
        # not implemented yet (multiple samples)
        pass

    def test_get_args(self):
        # not sure how to check arguments match inputs. Maybe can be left out
        # on the assumption that argparse is tested?
        pass

    def test_get_output_name(self):
        '''
        Checks that get_output_name() returns a correctly formatted filename
        string, either with the two sample names or with TPvsFP and the sample
        name if matching (i.e. happy & query of the same sample).
        '''
        self.assertEqual(
            FPTP.get_output_name(
                ['sample1-query.vcf', 'NA12878-query.vcf'], False
                ),
            "sample1_NA12878_QCdist.html"
            )
        self.assertEqual(
            FPTP.get_output_name(
                ['NA12878-test.vcf.gz', 'NA12878-query-test.vcf.gz']
                ),
            "TPvsFP_NA12878_QCdist.html"
            )

    def test_get_sample_names(self):
        '''
        Check that get_sample_names does not return an error
        '''
        self.assertEqual(
            FPTP.get_sample_names(
                'NA12878-test.vcf.gz',
                'NA12878-query-test.vcf.gz'
                ), None
            )
        self.assertRaises(
            AssertionError, FPTP.get_sample_names, 'NA12878-test.vcf.gz',
            'X219354_query-test.vcf.gz'
            )

    def test_infer_het_hom(self):
        '''
        Check that infer_het_hom() returns the appropriate string
        '''
        self.assertEqual(FPTP.infer_het_hom('1/1'), 'homalt')
        self.assertEqual(FPTP.infer_het_hom('0/1'), 'het')
        self.assertEqual(FPTP.infer_het_hom('0/0'), 'homref')

    def test_infer_snp_indel(self):
        '''
        Check that infer_snp_indel() returns the appropriate string
        '''
        self.assertEqual(FPTP.infer_snp_indel('AA', 'A'), 'INDEL')
        self.assertEqual(FPTP.infer_snp_indel('A', 'CAT'), 'INDEL')
        self.assertEqual(FPTP.infer_snp_indel('A', 'G'), 'SNP')

    def test_main(self):
        # do I need to test this?
        pass

    def test_make_lists(self):
        '''
        Check that make_lists() returns an appropriately filtered list for the
        requested metric.
        '''
        self.assertEqual(
            FPTP.make_lists(
                test_merged_dict, 'info_DP', ['TP', 'FP'],
                snp_indel='SNP', hethom='het'
            ), [['TP', 21.0, 43.0, 63.0], ['FP']]
            )

    def test_make_html(self):
        # apparently validating html is hard
        pass

    def test_make_plots(self):
        '''
        Check that make_plots() returns a plotly figure object matching the
        example
        '''
        # discussed with Matt - can't test something so big as no orthogonal
        # method for making the example
        pass

    def test_make_report(self):
        FPTP.make_report('Test_String', 'Test_file')
        with open('Test_file', 'r') as t:
            self.assertEqual(t.readlines()[0], 'Test_String')
        os.remove('Test_file')

    def test_make_tiled_figure(self):
        # TODO - need to make 4 input plots (or maybe just reuse the same
        # one 4 times...)
        pass

    def test_merge_happy_query(self):
        '''
        Check that merge_happy_query() merges two dictionaries correctly and
        returns a merged dict matching the example.
        '''
        self.assertEqual(
            FPTP.merge_happy_query(test_happy_dict, test_query_dict),
            test_merged_dict
            )

    def test_merge_samples(self):
        # not implemented yet (multiple samples)
        pass

    def test_parse_happy(self):
        '''
        Check that parse_happy() produces the correct dictionary structure from
        the test file (should match example).
        '''
        self.assertEqual(
            FPTP.parse_happy("NA12878-test.vcf.gz"), test_happy_dict
            )

    def test_parse_query(self):
        '''
        Check that parse_query() produces the correct dictionary structure from
        the test file (should match example).
        '''
        self.assertEqual(
            FPTP.parse_query("NA12878-query-test.vcf.gz"), test_query_dict
            )


if __name__ == '__main__':
    unittest.main()
