import FPTP
import unittest


1	949608	.	G	A,C	.	.	BS=949608;Regions=CONF,TS_contained	GT:BD:BK:QQ:BI:BVT:BLT	1|0:TP:gm:3222.77:ti:SNP:het	0/1:TP:gm:3222.77:ti:SNP:het
1	120612006	.	G	A	.	.	BS=120612006	GT:QQ:BD:BK:BI:BVT:BLT	.:.:UNK:.:.:NOCALL:nocall	0/1:2779.77:UNK:.:ti:SNP:het
2	238244863	.	TGCA	T	.	.	BS=238244864;Regions=CONF,TS_contained	GT:BD:BK:QQ:BI:BVT:BLT	1|0:TP:gm:3284.73:d1_5:INDEL:het	0/1:TP:gm:3284.73:d1_5:INDEL:het
6	161152240	.	G	A	.	.	BS=161152240	GT:QQ:BD:BK:BI:BVT:BLT	.:.:UNK:.:.:NOCALL:nocall	1/1:2892.77:UNK:.:ti:SNP:homalt
8	145738767	.	CG	C	.	.	BS=145738768	GT:QQ:BD:BK:BI:BVT:BLT	.:.:UNK:.:.:NOCALL:nocall	1/1:16296.7:UNK:.:d1_5:INDEL:homalt
11	68192690	.	G	A	.	.	BS=68192690;Regions=CONF,TS_contained	GT:BD:BK:QQ:BI:BVT:BLT	0|1:TP:gm:3303.77:ti:SNP:het	0/1:TP:gm:3303.77:ti:SNP:het
15	90210263	.	A	G	.	.	BS=90210263;Regions=CONF,TS_contained	GT:BD:BK:QQ:BI:BVT:BLT	0|1:TP:gm:2582.77:ti:SNP:het	0/1:TP:gm:2582.77:ti:SNP:het
16	89167443	.	T	C	.	.	BS=89167443;Regions=CONF,TS_contained	GT:BD:BK:QQ:BI:BVT:BLT	1|1:TP:gm:9255.77:ti:SNP:homalt	1/1:TP:gm:9255.77:ti:SNP:homalt
19	920642	.	T	A	.	.	BS=920642;Regions=CONF,TS_contained	GT:BD:BK:QQ:BI:BVT:BLT	1|1:TP:gm:8189.77:tv:SNP:homalt	1/1:TP:gm:8189.77:tv:SNP:homalt
21	47783796	.	T	C	.	.	BS=47783796;Regions=CONF,TS_contained	GT:BD:BK:QQ:BI:BVT:BLT	1|1:TP:gm:8294.77:ti:SNP:homalt	1/1:TP:gm:8294.77:ti:SNP:homalt

class TestModule(unittest.TestCase):
    def test_calculate_centiles(self):
        pass

    def test_check_happy_query_match(self):
        self.assertEqual(FPTP.check_happy_query_match("NA12878_test.vcf.gz", "NA12878-query-test.vcf.gz"),True)
        self.assertRaises(AssertionError, FPTP.check_happy_query_match, "NA12878_test.vcf.gz" "sample1-query.vcf")

    def test_check_metrics(self):
        pass

    def test_check_multiple_query_metrics(self):
        pass

    def test_create_plot(self):
        pass

    def test_decide_bins(self):
        pass

    def test_get_args(self):
        pass

    def test_get_output_name(self):
        pass

    def test_get_sample_names(self):
        pass

    def test_infer_het_hom(self):
        pass

    def test_infer_snp_indel(self):
        pass

    def test_main(self):
        pass

    def test_make_arrays(self):
        pass

    def test_make_html(self):
        pass

    def test_make_plots(self):
        pass

    def test_make_report(self):
        pass

    def test_make_subplots(self):
        pass

    def test_make_tiled_figure(self):
        pass

    def test_merge_happy_query(self):
        pass

    def test_merge_samples(self):
        pass

    def test_parse_happy(self):
        pass

    def test_parse_query(self):
        pass


if __name__ == '__main__':
    unittest.main()
