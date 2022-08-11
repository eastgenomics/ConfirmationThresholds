import FPTP
import unittest
from io import StringIO


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

test_fig = Figure({
    'data': [{'alignmentgroup': 'True',
              'bingroup': 'x',
              'hovertemplate': 'TPFP=TP<br>values=%{x}<br>count=%{y}<extra></extra>',
              'legendgroup': 'TP',
              'marker': {'color': '#636efa', 'opacity': 0.5, 'pattern': {'shape': ''}},
              'name': 'TP',
              'offsetgroup': 'TP',
              'orientation': 'v',
              'showlegend': True,
              'type': 'histogram',
              'x': array([0.749, 0.742, 0.772, ..., 0.758, 0.738, 0.864]),
              'xaxis': 'x',
              'yaxis': 'y'},
             {'alignmentgroup': 'True',
              'boxpoints': 'all',
              'customdata': array([[47.33904465],
                                   [45.01557632],
                                   [54.17964694],
                                   ...,
                                   [50.03894081],
                                   [43.83437175],
                                   [73.33852544]]),
              'fillcolor': 'rgba(255,255,255,0)',
              'hoveron': 'points',
              'hovertemplate': 'TPFP=TP<br>values=%{x}<br>centiles=%{customdata[0]}<extra></extra>',
              'jitter': 0,
              'legendgroup': 'TP',
              'line': {'color': 'rgba(255,255,255,0)'},
              'marker': {'color': '#636efa', 'symbol': 'line-ns-open'},
              'name': 'TP',
              'offsetgroup': 'TP',
              'showlegend': False,
              'type': 'box',
              'x': array([0.749, 0.742, 0.772, ..., 0.758, 0.738, 0.864]),
              'xaxis': 'x2',
              'yaxis': 'y2'},
             {'alignmentgroup': 'True',
              'bingroup': 'x',
              'hovertemplate': 'TPFP=FP<br>values=%{x}<br>count=%{y}<extra></extra>',
              'legendgroup': 'FP',
              'marker': {'color': '#EF553B', 'opacity': 0.5, 'pattern': {'shape': ''}},
              'name': 'FP',
              'offsetgroup': 'FP',
              'orientation': 'v',
              'showlegend': True,
              'type': 'histogram',
              'x': array([0.638, 0.792, 0.795, 0.824, 0.793, 0.832, 0.811, 0.772, 0.919, 1.02 ,
                          0.223, 0.788, 0.842, 0.876, 1.118, 5.025, 0.437, 1.125, 0.906, 1.131,
                          0.616, 0.663, 0.619, 1.185, 1.906, 0.582, 0.771, 1.849, 0.874, 0.56 ,
                          0.54 , 0.601, 0.424, 0.571, 0.671, 0.609, 1.18 , 0.787, 0.871, 0.809]),
              'xaxis': 'x',
              'yaxis': 'y'},
             {'alignmentgroup': 'True',
              'boxpoints': 'all',
              'customdata': array([[ 30. ],
                                   [ 47.5],
                                   [ 52.5],
                                   [ 60. ],
                                   [ 50. ],
                                   [ 62.5],
                                   [ 57.5],
                                   [ 40. ],
                                   [ 77.5],
                                   [ 80. ],
                                   [  2.5],
                                   [ 45. ],
                                   [ 65. ],
                                   [ 72.5],
                                   [ 82.5],
                                   [100. ],
                                   [  7.5],
                                   [ 85. ],
                                   [ 75. ],
                                   [ 87.5],
                                   [ 25. ],
                                   [ 32.5],
                                   [ 27.5],
                                   [ 92.5],
                                   [ 97.5],
                                   [ 17.5],
                                   [ 37.5],
                                   [ 95. ],
                                   [ 70. ],
                                   [ 12.5],
                                   [ 10. ],
                                   [ 20. ],
                                   [  5. ],
                                   [ 15. ],
                                   [ 35. ],
                                   [ 22.5],
                                   [ 90. ],
                                   [ 42.5],
                                   [ 67.5],
                                   [ 55. ]]),
              'fillcolor': 'rgba(255,255,255,0)',
              'hoveron': 'points',
              'hovertemplate': 'TPFP=FP<br>values=%{x}<br>centiles=%{customdata[0]}<extra></extra>',
              'jitter': 0,
              'legendgroup': 'FP',
              'line': {'color': 'rgba(255,255,255,0)'},
              'marker': {'color': '#EF553B', 'symbol': 'line-ns-open'},
              'name': 'FP',
              'offsetgroup': 'FP',
              'showlegend': False,
              'type': 'box',
              'x': array([0.638, 0.792, 0.795, 0.824, 0.793, 0.832, 0.811, 0.772, 0.919, 1.02 ,
                          0.223, 0.788, 0.842, 0.876, 1.118, 5.025, 0.437, 1.125, 0.906, 1.131,
                          0.616, 0.663, 0.619, 1.185, 1.906, 0.582, 0.771, 1.849, 0.874, 0.56 ,
                          0.54 , 0.601, 0.424, 0.571, 0.671, 0.609, 1.18 , 0.787, 0.871, 0.809]),
              'xaxis': 'x2',
              'yaxis': 'y2'}],
    'layout': {'barmode': 'overlay',
               'legend': {'title': {'text': 'TPFP'}, 'tracegroupgap': 0},
               'margin': {'t': 60},
               'template': '...',
               'xaxis': {'anchor': 'y', 'domain': [0.0, 1.0], 'title': {'text': 'values'}},
               'xaxis2': {'anchor': 'y2', 'domain': [0.0, 1.0], 'matches': 'x', 'showgrid': True, 'showticklabels': False},
               'yaxis': {'anchor': 'x', 'domain': [0.0, 0.7326], 'title': {'text': 'count'}},
               'yaxis2': {'anchor': 'x2',
                          'domain': [0.7426, 1.0],
                          'matches': 'y2',
                          'showgrid': False,
                          'showline': False,
                          'showticklabels': False,
                          'ticks': ''}}
})



class TestModule(unittest.TestCase):
    def test_calculate_centiles(self):
        pass

    def test_check_happy_query_match(self):
        self.assertEqual(FPTP.check_happy_query_match("NA12878_test.vcf.gz", "NA12878-query-test.vcf.gz"),True)
        self.assertRaises(AssertionError, FPTP.check_happy_query_match, "NA12878_test.vcf.gz" "sample1-query.vcf")

    def test_check_metrics(self):
        self.assertEqual(FPTP.check_metrics("NA12878-query-test.vcf.gz",['AD','DP','GQ','AC','AF','AN','BaseQRankSum']),[['AD','DP','GQ']['AC','AF','AN','BaseQRankSum']])
        self.assertEqual(FPTP.check_metrics("NA12878-query-test.vcf.gz",['AD','DP','GQ','AC','AF','AN','BaseQRankSum','MQ']),[['AD','DP','GQ']['AC','AF','AN','BaseQRankSum']])

    def test_check_multiple_query_metrics(self):
        pass

    def test_create_plot(self):
        pass

    def test_decide_bins(self):
        pass

    def test_get_args(self):
        pass

    def test_get_output_name(self):
        self.assertEqual(FPTP.get_output_name(['sample1-query.vcf', 'NA12878-query.vcf'], False),"sample1_NA12878_QCDist.html")
        self.assertEqual(FPTP.get_output_name(['NA12878-query-test.vcf.gz']),"TPvsFP_NA12878_QCDist.html")

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
        self.assertEqual(FPTP.merge_samples(test_dict_1, test_dict_2), test_dict_3)

    def test_parse_happy(self):
        test_dict = {}
        test_dict['1_949608_A_C'] = ['TP','SNP','het']
        test_dict['1_949654_A_G'] = ['TP','SNP','homalt']
        test_dict['1_981931_A_G'] = ['TP','SNP','het']
        test_dict['1_982994_T_C'] = ['TP','SNP','het']
        self.assertEqual(FPTP.parse_happy("NA12878_test.vcf.gz"), test_dict)

    def test_parse_query(self):
        test_dict = {}
        test_dict['1_229673_A_C'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'3', 'format_GT':'1/1', 'format_AD':'0,3', 'format_DP':'3'}
        test_dict['1_778302_C_CCT'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict['1_787262_C_G'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'2', 'format_GT':'1/1', 'format_AD':'0,2', 'format_DP':'2'}
        test_dict['1_787399_G_T'] = {'info_AC':'2', 'info_AF':'1', 'info_AN':'2', 'info_DB':True, 'info_DP':'10', 'format_GT':'1/1', 'format_AD':'0,10', 'format_DP':'10'}
        self.assertEqual(FPTP.parse_query("NA12878-query-test.vcf.gz"), test_dict)


if __name__ == '__main__':
    unittest.main()