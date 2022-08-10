from FPTP import *
import unittest
from io import StringIO


class TestRevComp(unittest.TestCase):
    def test_parseArgs(self):
        # not sure how to test this one
        pass


    def test_check_happy_query_match(self):
        self.assertEqual(check_happy_query_match("NA12878-happy.vcf", "NA12878-query.vcf"),True)
        self.assertRaises(AssertionError, check_happy_query_match, "NA12878-happy.vcf" "sample1-query.vcf")


    def test_check_multiple_query_metrics(self):
        self.assertEqual(check_multiple_query_metrics(["NA12878-query.vcf", "sample1-query.vcf"],['AD','DP','GQ','AC','AF','AN','BaseQRankSum']),[['AD','DP','GQ']['AC','AF','AN','BaseQRankSum']])
        self.assertEqual(check_multiple_query_metrics(["NA12878-query.vcf", "sample1-query.vcf"],['AD','DP','GQ','AC','AF','AN','BaseQRankSum','MQ']),[['AD','DP','GQ']['AC','AF','AN','BaseQRankSum']])


    def test_check_metrics(self):
        self.assertEqual(check_metrics("sample1-query.vcf",['AD','DP','GQ','AC','AF','AN','BaseQRankSum']),[['AD','DP','GQ']['AC','AF','AN','BaseQRankSum']])
        self.assertEqual(check_metrics("sample1-query.vcf",['AD','DP','GQ','AC','AF','AN','BaseQRankSum','MQ']),[['AD','DP','GQ']['AC','AF','AN','BaseQRankSum']])


    def test_parse_query(self):
        test_dict = {}
        test_dict['1_229673_A_C'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'3', 'format_GT':'1/1', 'format_AD':'0,3', 'format_DP':'3'}
        test_dict['1_778302_C_CCT'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict['1_787262_C_G'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict['1_787399_G_T'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'10', 'format_GT':'1/1', 'format_AD':'0,10', 'format_DP':'10'}
        self.assertEqual(parse_query("sample1-query-short.vcf"), test_dict)


    def test_parse_happy(self):
        test_dict = {}
        test_dict['1_949608_A_C'] = ['TP','SNP','het']
        test_dict['1_949654_A_G'] = ['TP','SNP','homalt']
        test_dict['1_981931_A_G'] = ['TP','SNP','het']
        test_dict['1_982994_T_C'] = ['TP','SNP','het']
        self.assertEqual(parse_happy("NA12878-happy-short.vcf"), test_dict)


    def test_get_output_name(self):
        self.assertEqual(get_output_name(['sample1-query.vcf', 'NA12878-query.vcf'], False),"sample1_NA12878_QCDist.html")
        self.assertEqual(get_output_name(['NA12878-query.vcf']),"TPvsFP_NA12878_QCDist.html")


    def test_makeReport(self):
        pass


    def test_createPlot(self):
        pass


    def test_mergeSamples(self):
        test_dict_1 = {}
        test_dict_1['1_229673_A_C'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'3', 'format_GT':'1/1', 'format_AD':'0,3', 'format_DP':'3'}
        test_dict_1['1_778302_C_CCT'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict_1['1_787262_C_G'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict_1['1_787399_G_T'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'10', 'format_GT':'1/1', 'format_AD':'0,10', 'format_DP':'10'}
        test_dict_2 = {}
        test_dict_2['1_329673_G_C'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'3', 'format_GT':'1/1', 'format_AD':'0,3', 'format_DP':'3'}
        test_dict_2['1_478302_C_CT'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict_2['1_587262_C_T'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict_2['1_687399_A_G'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'10', 'format_GT':'1/1', 'format_AD':'0,10', 'format_DP':'10'}
        test_dict_3 = {}
        test_dict_3['1_229673_A_C'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'3', 'format_GT':'1/1', 'format_AD':'0,3', 'format_DP':'3'}
        test_dict_3['1_778302_C_CCT'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict_3['1_787262_C_G'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict_3['1_787399_G_T'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'10', 'format_GT':'1/1', 'format_AD':'0,10', 'format_DP':'10'}
        test_dict_3['1_329673_G_C'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'3', 'format_GT':'1/1', 'format_AD':'0,3', 'format_DP':'3'}
        test_dict_3['1_478302_C_CT'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict_3['1_587262_C_T'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict_3['1_687399_A_G'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'10', 'format_GT':'1/1', 'format_AD':'0,10', 'format_DP':'10'}
        self.assertEqual(mergeSamples(test_dict_1, test_dict_2), test_dict_3)


    def test_makeArrays(self):
        pass


if __name__ == '__main__':
    unittest.main()