from FPTP import *
import unittest
from io import StringIO


class TestRevComp(unittest.TestCase):
    def test_parseArgs(self):
        # not sure how to test this one
        pass


    def test_checkHappyQueryMatch(self):
        self.assertEqual(checkHappyQueryMatch("NA12878-happy.vcf", "NA12878-query.vcf"),True)
        self.assertRaises(AssertionError, checkHappyQueryMatch, "NA12878-happy.vcf" "sample1-query.vcf")


    def test_checkMultipleQueryMetrics(self):
        self.assertEqual(checkMultipleQueryMetrics(["NA12878-query.vcf", "sample1-query.vcf"],['AD','DP','GQ','AC','AF','AN','BaseQRankSum']),[['AD','DP','GQ']['AC','AF','AN','BaseQRankSum']])
        self.assertEqual(checkMultipleQueryMetrics(["NA12878-query.vcf", "sample1-query.vcf"],['AD','DP','GQ','AC','AF','AN','BaseQRankSum','MQ']),[['AD','DP','GQ']['AC','AF','AN','BaseQRankSum']])


    def test_checkMetrics(self):
        self.assertEqual(checkMetrics("sample1-query.vcf",['AD','DP','GQ','AC','AF','AN','BaseQRankSum']),[['AD','DP','GQ']['AC','AF','AN','BaseQRankSum']])
        self.assertEqual(checkMetrics("sample1-query.vcf",['AD','DP','GQ','AC','AF','AN','BaseQRankSum','MQ']),[['AD','DP','GQ']['AC','AF','AN','BaseQRankSum']])


    def test_parseQuery(self):
        test_dict = {}
        test_dict['1_229673_A_C'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'3', 'format_GT':'1/1', 'format_AD':'0,3', 'format_DP':'3'}
        test_dict['1_778302_C_CCT'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict['1_787262_C_G'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict['1_787399_G_T'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'10', 'format_GT':'1/1', 'format_AD':'0,10', 'format_DP':'10'}
        self.assertEqual(parseQuery("sample1-query-short.vcf"), test_dict)


    def test_parseHappy(self):
        test_dict = {}
        test_dict['1_949608_A_C'] = ['TP','SNP','het']
        test_dict['1_949654_A_G'] = ['TP','SNP','homalt']
        test_dict['1_981931_A_G'] = ['TP','SNP','het']
        test_dict['1_982994_T_C'] = ['TP','SNP','het']
        self.assertEqual(parseHappy("NA12878-happy-short.vcf"), test_dict)


    def test_getOutputName(self):
        self.assertEqual(getOutputName(['sample1-query.vcf', 'NA12878-query.vcf'], False),"sample1_NA12878_QCDist.html")
        self.assertEqual(getOutputName(['NA12878-query.vcf']),"TPvsFP_NA12878_QCDist.html")


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