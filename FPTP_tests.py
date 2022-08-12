from cgi import test
import FPTP
import unittest
from io import StringIO

# these test objects are correct as long as the test input files bundled with the module are unchanged
test_query_dict = {'1_949608_G_A': {'info_AC': 2, 'info_AF': 1.0, 'info_AN': 2, 'info_DP': 21, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 2, 'info_MLEAF': 1.0, 'info_MQ': 26.41, 'info_QD': 28.85, 'info_SOR': 6.184, 'format_GT': '1/1', 'format_AD': [0, 21], 'format_DP': 21, 'format_GQ': 63, 'format_PL': [634, 63, 0]}, '1_120612006_G_A': {'info_AC': 2, 'info_AF': 1.0, 'info_AN': 2, 'info_BaseQRankSum': 1.022, 'info_ClippingRankSum': 0.0, 'info_DB': True, 'info_DP': 84, 'info_ExcessHet': 3.0103, 'info_FS': 6.938, 'info_MLEAC': 2, 'info_MLEAF': 1.0, 'info_MQ': 25.88, 'info_MQRankSum': -1.034, 'info_QD': 42.34, 'info_ReadPosRankSum': -0.664, 'info_SOR': 2.788, 'format_GT': '1/1', 'format_AD': [1, 83], 'format_DP': 84, 'format_GQ': 99, 'format_PGT': '1|1', 'format_PID': '13116_T_G', 'format_PL': [3585, 189, 0]}, '1_14574_A_G': {'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': -4.914, 'info_ClippingRankSum': 0.0, 'info_DP': 61, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 27.88, 'info_MQRankSum': -6.387, 'info_QD': 7.63, 'info_ReadPosRankSum': -0.422, 'info_SOR': 0.298, 'format_GT': '0/1', 'format_AD': [38, 22], 'format_DP': 60, 'format_GQ': 99, 'format_PL': [486, 0, 1058]}, '1_13838_C_T': {'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': -0.687, 'info_ClippingRankSum': 0.0, 'info_DB': True, 'info_DP': 35, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 28.39, 'info_MQRankSum': -1.352, 'info_QD': 2.36, 'info_ReadPosRankSum': 0.049, 'info_SOR': 0.148, 'format_GT': '0/1', 'format_AD': [30, 5], 'format_DP': 35, 'format_GQ': 99, 'format_PGT': '0|1', 'format_PID': '13813_T_G', 'format_PL': [111, 0, 1412]}, '2_238244863_TGCA_T': {'info_AC': 2, 'info_AF': 1.0, 'info_AN': 2, 'info_BaseQRankSum': -0.752, 'info_ClippingRankSum': 0.0, 'info_DB': True, 'info_DP': 79, 'info_ExcessHet': 3.0103, 'info_FS': 6.935, 'info_MLEAC': 2, 'info_MLEAF': 1.0, 'info_MQ': 25.57, 'info_MQRankSum': -1.147, 'info_QD': 45.02, 'info_ReadPosRankSum': -0.62, 'info_SOR': 2.788, 'format_GT': '1/1', 'format_AD': [1, 78], 'format_DP': 79, 'format_GQ': 99, 'format_PGT': '1|1', 'format_PID': '13116_T_G', 'format_PL': [3585, 189, 0]}, '6_161152240_G_A': {'info_AC': 2, 'info_AF': 1.0, 'info_AN': 2, 'info_BaseQRankSum': -3.811, 'info_ClippingRankSum': 0.0, 'info_DP': 107, 'info_ExcessHet': 3.0103, 'info_FS': 43.338, 'info_MLEAC': 2, 'info_MLEAF': 1.0, 'info_MQ': 23.82, 'info_MQRankSum': 0.165, 'info_QD': 38.17, 'info_ReadPosRankSum': -0.629, 'info_SOR': 3.196, 'format_GT': '1/1', 'format_AD': [6, 99], 'format_DP': 105, 'format_GQ': 42, 'format_PL': [4045, 42, 0]}, '8_145738767_CG_C': {'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': 1.604, 'info_ClippingRankSum': 0.0, 'info_DP': 43, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 27.39, 'info_MQRankSum': -0.791, 'info_QD': 1.51, 'info_ReadPosRankSum': 1.061, 'info_SOR': 0.095, 'format_GT': '0/1', 'format_AD': [38, 5], 'format_DP': 43, 'format_GQ': 93, 'format_PGT': '0|1', 'format_PID': '13813_T_G', 'format_PL': [93, 0, 1577]}, '11_68192690_G_A': {'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': 3.147, 'info_ClippingRankSum': 0.0, 'info_DP': 43, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 28.36, 'info_MQRankSum': -5.35, 'info_QD': 7.34, 'info_ReadPosRankSum': -0.285, 'info_SOR': 0.223, 'format_GT': '0/1', 'format_AD': [29, 14], 'format_DP': 43, 'format_GQ': 99, 'format_PL': [344, 0, 841]}, '15_90210263_A_G': {'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': 4.908, 'info_ClippingRankSum': 0.0, 'info_DP': 63, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 27.77, 'info_MQRankSum': -3.548, 'info_QD': 0.55, 'info_ReadPosRankSum': -0.167, 'info_SOR': 0.033, 'format_GT': '0/1', 'format_AD': [54, 9], 'format_DP': 63, 'format_GQ': 63, 'format_PL': [63, 0, 2016]}, '16_89167443_T_C': {'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': -5.087, 'info_ClippingRankSum': 0.0, 'info_DP': 64, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 27.76, 'info_MQRankSum': -6.212, 'info_QD': 11.97, 'info_ReadPosRankSum': 0.467, 'info_SOR': 0.252, 'format_GT': '0/1', 'format_AD': [42, 22], 'format_DP': 64, 'format_GQ': 99, 'format_PGT': '0|1', 'format_PID': '14599_T_A', 'format_PL': [794, 0, 1739]}, '19_920642_T_A': {'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': 5.779, 'info_ClippingRankSum': 0.0, 'info_DP': 66, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 27.74, 'info_MQRankSum': -6.312, 'info_QD': 13.86, 'info_ReadPosRankSum': 0.783, 'info_SOR': 0.232, 'format_GT': '0/1', 'format_AD': [44, 22], 'format_DP': 66, 'format_GQ': 99, 'format_PGT': '0|1', 'format_PID': '14599_T_A', 'format_PL': [943, 0, 1769]}, '21_47783796_T_C': {'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': 5.661, 'info_ClippingRankSum': 0.0, 'info_DP': 75, 'info_ExcessHet': 3.0103, 'info_FS': 4.127, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 28.86, 'info_MQRankSum': -4.114, 'info_QD': 13.37, 'info_ReadPosRankSum': -0.719, 'info_SOR': 0.091, 'format_GT': '0/1', 'format_AD': [46, 29], 'format_DP': 75, 'format_GQ': 99, 'format_PGT': '0|1', 'format_PID': '14599_T_A', 'format_PL': [1031, 0, 1816]}}
test_happy_dict = {'1_949608_G_A': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'het'}, '1_120612006_G_A': {'TPFP_or_samplename': 'UNK', 'snp_indel': 'SNP', 'HETHOM': 'het'}, '2_238244863_TGCA_T': {'TPFP_or_samplename': 'TP', 'snp_indel': 'INDEL', 'HETHOM': 'het'}, '6_161152240_G_A': {'TPFP_or_samplename': 'UNK', 'snp_indel': 'SNP', 'HETHOM': 'homalt'}, '8_145738767_CG_C': {'TPFP_or_samplename': 'UNK', 'snp_indel': 'INDEL', 'HETHOM': 'homalt'}, '11_68192690_G_A': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'het'}, '15_90210263_A_G': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'het'}, '16_89167443_T_C': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'homalt'}, '19_920642_T_A': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'homalt'}, '21_47783796_T_C': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'homalt'}}
test_merged_dict = {'1_949608_G_A': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'het', 'info_AC': 2, 'info_AF': 1.0, 'info_AN': 2, 'info_DP': 21, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 2, 'info_MLEAF': 1.0, 'info_MQ': 26.41, 'info_QD': 28.85, 'info_SOR': 6.184, 'format_GT': '1/1', 'format_AD': [0, 21], 'format_DP': 21, 'format_GQ': 63, 'format_PL': [634, 63, 0]}, '1_120612006_G_A': {'TPFP_or_samplename': 'UNK', 'snp_indel': 'SNP', 'HETHOM': 'het', 'info_AC': 2, 'info_AF': 1.0, 'info_AN': 2, 'info_BaseQRankSum': 1.022, 'info_ClippingRankSum': 0.0, 'info_DB': True, 'info_DP': 84, 'info_ExcessHet': 3.0103, 'info_FS': 6.938, 'info_MLEAC': 2, 'info_MLEAF': 1.0, 'info_MQ': 25.88, 'info_MQRankSum': -1.034, 'info_QD': 42.34, 'info_ReadPosRankSum': -0.664, 'info_SOR': 2.788, 'format_GT': '1/1', 'format_AD': [1, 83], 'format_DP': 84, 'format_GQ': 99, 'format_PGT': '1|1', 'format_PID': '13116_T_G', 'format_PL': [3585, 189, 0]}, '2_238244863_TGCA_T': {'TPFP_or_samplename': 'TP', 'snp_indel': 'INDEL', 'HETHOM': 'het', 'info_AC': 2, 'info_AF': 1.0, 'info_AN': 2, 'info_BaseQRankSum': -0.752, 'info_ClippingRankSum': 0.0, 'info_DB': True, 'info_DP': 79, 'info_ExcessHet': 3.0103, 'info_FS': 6.935, 'info_MLEAC': 2, 'info_MLEAF': 1.0, 'info_MQ': 25.57, 'info_MQRankSum': -1.147, 'info_QD': 45.02, 'info_ReadPosRankSum': -0.62, 'info_SOR': 2.788, 'format_GT': '1/1', 'format_AD': [1, 78], 'format_DP': 79, 'format_GQ': 99, 'format_PGT': '1|1', 'format_PID': '13116_T_G', 'format_PL': [3585, 189, 0]}, '6_161152240_G_A': {'TPFP_or_samplename': 'UNK', 'snp_indel': 'SNP', 'HETHOM': 'homalt', 'info_AC': 2, 'info_AF': 1.0, 'info_AN': 2, 'info_BaseQRankSum': -3.811, 'info_ClippingRankSum': 0.0, 'info_DP': 107, 'info_ExcessHet': 3.0103, 'info_FS': 43.338, 'info_MLEAC': 2, 'info_MLEAF': 1.0, 'info_MQ': 23.82, 'info_MQRankSum': 0.165, 'info_QD': 38.17, 'info_ReadPosRankSum': -0.629, 'info_SOR': 3.196, 'format_GT': '1/1', 'format_AD': [6, 99], 'format_DP': 105, 'format_GQ': 42, 'format_PL': [4045, 42, 0]}, '8_145738767_CG_C': {'TPFP_or_samplename': 'UNK', 'snp_indel': 'INDEL', 'HETHOM': 'homalt', 'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': 1.604, 'info_ClippingRankSum': 0.0, 'info_DP': 43, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 27.39, 'info_MQRankSum': -0.791, 'info_QD': 1.51, 'info_ReadPosRankSum': 1.061, 'info_SOR': 0.095, 'format_GT': '0/1', 'format_AD': [38, 5], 'format_DP': 43, 'format_GQ': 93, 'format_PGT': '0|1', 'format_PID': '13813_T_G', 'format_PL': [93, 0, 1577]}, '11_68192690_G_A': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'het', 'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': 3.147, 'info_ClippingRankSum': 0.0, 'info_DP': 43, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 28.36, 'info_MQRankSum': -5.35, 'info_QD': 7.34, 'info_ReadPosRankSum': -0.285, 'info_SOR': 0.223, 'format_GT': '0/1', 'format_AD': [29, 14], 'format_DP': 43, 'format_GQ': 99, 'format_PL': [344, 0, 841]}, '15_90210263_A_G': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'het', 'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': 4.908, 'info_ClippingRankSum': 0.0, 'info_DP': 63, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 27.77, 'info_MQRankSum': -3.548, 'info_QD': 0.55, 'info_ReadPosRankSum': -0.167, 'info_SOR': 0.033, 'format_GT': '0/1', 'format_AD': [54, 9], 'format_DP': 63, 'format_GQ': 63, 'format_PL': [63, 0, 2016]}, '16_89167443_T_C': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'homalt', 'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': -5.087, 'info_ClippingRankSum': 0.0, 'info_DP': 64, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 27.76, 'info_MQRankSum': -6.212, 'info_QD': 11.97, 'info_ReadPosRankSum': 0.467, 'info_SOR': 0.252, 'format_GT': '0/1', 'format_AD': [42, 22], 'format_DP': 64, 'format_GQ': 99, 'format_PGT': '0|1', 'format_PID': '14599_T_A', 'format_PL': [794, 0, 1739]}, '19_920642_T_A': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'homalt', 'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': 5.779, 'info_ClippingRankSum': 0.0, 'info_DP': 66, 'info_ExcessHet': 3.0103, 'info_FS': 0.0, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 27.74, 'info_MQRankSum': -6.312, 'info_QD': 13.86, 'info_ReadPosRankSum': 0.783, 'info_SOR': 0.232, 'format_GT': '0/1', 'format_AD': [44, 22], 'format_DP': 66, 'format_GQ': 99, 'format_PGT': '0|1', 'format_PID': '14599_T_A', 'format_PL': [943, 0, 1769]}, '21_47783796_T_C': {'TPFP_or_samplename': 'TP', 'snp_indel': 'SNP', 'HETHOM': 'homalt', 'info_AC': 1, 'info_AF': 0.5, 'info_AN': 2, 'info_BaseQRankSum': 5.661, 'info_ClippingRankSum': 0.0, 'info_DP': 75, 'info_ExcessHet': 3.0103, 'info_FS': 4.127, 'info_MLEAC': 1, 'info_MLEAF': 0.5, 'info_MQ': 28.86, 'info_MQRankSum': -4.114, 'info_QD': 13.37, 'info_ReadPosRankSum': -0.719, 'info_SOR': 0.091, 'format_GT': '0/1', 'format_AD': [46, 29], 'format_DP': 75, 'format_GQ': 99, 'format_PGT': '0|1', 'format_PID': '14599_T_A', 'format_PL': [1031, 0, 1816]}}
test_arrays = [['TP', 21.0, 79.0, 43.0, 63.0, 64.0, 66.0, 75.0], ['FP']]
test_figure = ('Figure({\
    \'data\': [{\'alignmentgroup\': \'True\',\
              \'bingroup\': \'x\',\
              \'hovertemplate\': \'True/False Positive=TP<br>Bin=%{x}<extra></extra>\',\
              \'legendgroup\': \'TP\',\
              \'marker\': {\'color\': \'#636efa\', \'opacity\': 0.5, \'pattern\': {\'shape\': \'\'}},\
              \'name\': \'TP\',\
              \'offsetgroup\': \'TP\',\
              \'orientation\': \'v\',\
              \'showlegend\': True,\
              \'type\': \'histogram\',\
              \'x\': array([21., 79., 43., 63., 64., 66., 75.]),\
              \'xaxis\': \'x\',\
              \'yaxis\': \'y\'},\
             {\'alignmentgroup\': \'True\',\
              \'boxpoints\': \'all\',\
              \'customdata\': array([[ 14.29],\
                                   [100.  ],\
                                   [ 28.57],\
                                   [ 42.86],\
                                   [ 57.14],\
                                   [ 71.43],\
                                   [ 85.71]]),\
              \'fillcolor\': \'rgba(255,255,255,0)\',\
              \'hoveron\': \'points\',\
              \'hovertemplate\': \'<br>Metric value=%{x}<br>Centile=%{customdata[0]}<br><extra></extra>\',\
              \'jitter\': 0,\
              \'legendgroup\': \'TP\',\
              \'line\': {\'color\': \'rgba(255,255,255,0)\'},\
              \'marker\': {\'color\': \'#636efa\', \'symbol\': \'line-ns-open\'},\
              \'name\': \'TP\',\
              \'offsetgroup\': \'TP\',\
              \'showlegend\': False,\
              \'type\': \'box\',\
              \'x\': array([21., 79., 43., 63., 64., 66., 75.]),\
              \'xaxis\': \'x2\',\
              \'yaxis\': \'y2\'}],\
    \'layout\': {\'barmode\': \'overlay\',\
               \'legend\': {\'title\': {\'text\': \'TPFP\'}, \'tracegroupgap\': 0},\
               \'margin\': {\'t\': 60},\
               \'template\': \'...\',\
               \'xaxis\': {\'anchor\': \'y\', \'domain\': [0.0, 1.0], \'title\': {\'text\': \'values\'}},\
               \'xaxis2\': {\'anchor\': \'y2\', \'domain\': [0.0, 1.0], \'matches\': \'x\', \'showgrid\': True, \'showticklabels\': False},\
               \'yaxis\': {\'anchor\': \'x\', \'domain\': [0.0, 0.7326], \'title\': {\'text\': \'count\'}},\
               \'yaxis2\': {\'anchor\': \'x2\',\
                          \'domain\': [0.7426, 1.0],\
                          \'matches\': \'y2\',\
                          \'showgrid\': False,\
                          \'showline\': False,\
                          \'showticklabels\': False,\
                          \'ticks\': \'\'}}\
})\
')
test_tiled_figure = ('Figure({\
    \'data\': [{\'alignmentgroup\': \'True\',\
              \'bingroup\': \'x\',\
              \'hovertemplate\': \'True/False Positive=TP<br>Bin=%{x}<extra></extra>\',\
              \'legendgroup\': \'TP\',\
              \'marker\': {\'color\': \'#636efa\', \'opacity\': 0.5, \'pattern\': {\'shape\': \'\'}},\
              \'name\': \'TP\',\
              \'offsetgroup\': \'TP\',\
              \'orientation\': \'v\',\
              \'showlegend\': True,\
              \'type\': \'histogram\',\
              \'x\': array([2., 2., 2., 2., 2., 2., 2.]),\
              \'xaxis\': \'x\',\
              \'yaxis\': \'y\'},\
             {\'alignmentgroup\': \'True\',\
              \'boxpoints\': \'all\',\
              \'customdata\': array([[57.14],\
                                   [57.14],\
                                   [57.14],\
                                   [57.14],\
                                   [57.14],\
                                   [57.14],\
                                   [57.14]]),\
              \'fillcolor\': \'rgba(255,255,255,0)\',\
              \'hoveron\': \'points\',\
              \'hovertemplate\': \'<br>Metric value=%{x}<br>Centile=%{customdata[0]}<br><extra></extra>\',\
              \'jitter\': 0,\
              \'legendgroup\': \'TP\',\
              \'line\': {\'color\': \'rgba(255,255,255,0)\'},\
              \'marker\': {\'color\': \'#636efa\', \'symbol\': \'line-ns-open\'},\
              \'name\': \'TP\',\
              \'offsetgroup\': \'TP\',\
              \'showlegend\': False,\
              \'type\': \'box\',\
              \'x\': array([2., 2., 2., 2., 2., 2., 2.]),\
              \'xaxis\': \'x\',\
              \'yaxis\': \'y\'},\
             {\'alignmentgroup\': \'True\',\
              \'bingroup\': \'x\',\
              \'hovertemplate\': \'True/False Positive=TP<br>Bin=%{x}<extra></extra>\',\
              \'legendgroup\': \'TP\',\
              \'marker\': {\'color\': \'#636efa\', \'opacity\': 0.5, \'pattern\': {\'shape\': \'\'}},\
              \'name\': \'TP\',\
              \'offsetgroup\': \'TP\',\
              \'orientation\': \'v\',\
              \'showlegend\': True,\
              \'type\': \'histogram\',\
              \'x\': array([2., 2., 2., 2., 2., 2., 2.]),\
              \'xaxis\': \'x2\',\
              \'yaxis\': \'y2\'},\
             {\'alignmentgroup\': \'True\',\
              \'boxpoints\': \'all\',\
              \'customdata\': array([[57.14],\
                                   [57.14],\
                                   [57.14],\
                                   [57.14],\
                                   [57.14],\
                                   [57.14],\
                                   [57.14]]),\
              \'fillcolor\': \'rgba(255,255,255,0)\',\
              \'hoveron\': \'points\',\
              \'hovertemplate\': \'<br>Metric value=%{x}<br>Centile=%{customdata[0]}<br><extra></extra>\',\
              \'jitter\': 0,\
              \'legendgroup\': \'TP\',\
              \'line\': {\'color\': \'rgba(255,255,255,0)\'},\
              \'marker\': {\'color\': \'#636efa\', \'symbol\': \'line-ns-open\'},\
              \'name\': \'TP\',\
              \'offsetgroup\': \'TP\',\
              \'showlegend\': False,\
              \'type\': \'box\',\
              \'x\': array([2., 2., 2., 2., 2., 2., 2.]),\
              \'xaxis\': \'x2\',\
              \'yaxis\': \'y2\'},\
             {\'alignmentgroup\': \'True\',\
              \'bingroup\': \'x\',\
              \'hovertemplate\': \'True/False Positive=TP<br>Bin=%{x}<extra></extra>\',\
              \'legendgroup\': \'TP\',\
              \'marker\': {\'color\': \'#636efa\', \'opacity\': 0.5, \'pattern\': {\'shape\': \'\'}},\
              \'name\': \'TP\',\
              \'offsetgroup\': \'TP\',\
              \'orientation\': \'v\',\
              \'showlegend\': True,\
              \'type\': \'histogram\',\
              \'x\': array([2., 2., 2., 2.]),\
              \'xaxis\': \'x3\',\
              \'yaxis\': \'y3\'},\
             {\'alignmentgroup\': \'True\',\
              \'boxpoints\': \'all\',\
              \'customdata\': array([[62.5],\
                                   [62.5],\
                                   [62.5],\
                                   [62.5]]),\
              \'fillcolor\': \'rgba(255,255,255,0)\',\
              \'hoveron\': \'points\',\
              \'hovertemplate\': \'<br>Metric value=%{x}<br>Centile=%{customdata[0]}<br><extra></extra>\',\
              \'jitter\': 0,\
              \'legendgroup\': \'TP\',\
              \'line\': {\'color\': \'rgba(255,255,255,0)\'},\
              \'marker\': {\'color\': \'#636efa\', \'symbol\': \'line-ns-open\'},\
              \'name\': \'TP\',\
              \'offsetgroup\': \'TP\',\
              \'showlegend\': False,\
              \'type\': \'box\',\
              \'x\': array([2., 2., 2., 2.]),\
              \'xaxis\': \'x3\',\
              \'yaxis\': \'y3\'},\
             {\'alignmentgroup\': \'True\',\
              \'bingroup\': \'x\',\
              \'hovertemplate\': \'True/False Positive=TP<br>Bin=%{x}<extra></extra>\',\
              \'legendgroup\': \'TP\',\
              \'marker\': {\'color\': \'#636efa\', \'opacity\': 0.5, \'pattern\': {\'shape\': \'\'}},\
              \'name\': \'TP\',\
              \'offsetgroup\': \'TP\',\
              \'orientation\': \'v\',\
              \'showlegend\': True,\
              \'type\': \'histogram\',\
              \'x\': array([2., 2., 2.]),\
              \'xaxis\': \'x4\',\
              \'yaxis\': \'y4\'},\
             {\'alignmentgroup\': \'True\',\
              \'boxpoints\': \'all\',\
              \'customdata\': array([[66.67],\
                                   [66.67],\
                                   [66.67]]),\
              \'fillcolor\': \'rgba(255,255,255,0)\',\
              \'hoveron\': \'points\',\
              \'hovertemplate\': \'<br>Metric value=%{x}<br>Centile=%{customdata[0]}<br><extra></extra>\',\
              \'jitter\': 0,\
              \'legendgroup\': \'TP\',\
              \'line\': {\'color\': \'rgba(255,255,255,0)\'},\
              \'marker\': {\'color\': \'#636efa\', \'symbol\': \'line-ns-open\'},\
              \'name\': \'TP\',\
              \'offsetgroup\': \'TP\',\
              \'showlegend\': False,\
              \'type\': \'box\',\
              \'x\': array([2., 2., 2.]),\
              \'xaxis\': \'x4\',\
              \'yaxis\': \'y4\'}],\
    \'layout\': {\'annotations\': [{\'font\': {\'size\': 16},\
                                \'showarrow\': False,\
                                \'text\': \'SNP\',\
                                \'x\': 0.10625,\
                                \'xanchor\': \'center\',\
                                \'xref\': \'paper\',\
                                \'y\': 1.0,\
                                \'yanchor\': \'bottom\',\
                                \'yref\': \'paper\'},\
                               {\'font\': {\'size\': 16},\
                                \'showarrow\': False,\
                                \'text\': \'INDEL\',\
                                \'x\': 0.36875,\
                                \'xanchor\': \'center\',\
                                \'xref\': \'paper\',\
                                \'y\': 1.0,\
                                \'yanchor\': \'bottom\',\
                                \'yref\': \'paper\'},\
                               {\'font\': {\'size\': 16},\
                                \'showarrow\': False,\
                                \'text\': \'HET\',\
                                \'x\': 0.6312500000000001,\
                                \'xanchor\': \'center\',\
                                \'xref\': \'paper\',\
                                \'y\': 1.0,\
                                \'yanchor\': \'bottom\',\
                                \'yref\': \'paper\'},\
                               {\'font\': {\'size\': 16},\
                                \'showarrow\': False,\
                                \'text\': \'HOM\',\
                                \'x\': 0.89375,\
                                \'xanchor\': \'center\',\
                                \'xref\': \'paper\',\
                                \'y\': 1.0,\
                                \'yanchor\': \'bottom\',\
                                \'yref\': \'paper\'}],\
               \'height\': 500,\
               \'hovermode\': \'x unified\',\
               \'showlegend\': False,\
               \'template\': \'...\',\
               \'title\': {\'text\': \'info_AN\'},\
               \'width\': 1800,\
               \'xaxis\': {\'anchor\': \'y\', \'domain\': [0.0, 0.2125]},\
               \'xaxis2\': {\'anchor\': \'y2\', \'domain\': [0.2625, 0.475]},\
               \'xaxis3\': {\'anchor\': \'y3\', \'domain\': [0.525, 0.7375]},\
               \'xaxis4\': {\'anchor\': \'y4\', \'domain\': [0.7875, 1.0]},\
               \'yaxis\': {\'anchor\': \'x\', \'domain\': [0.0, 1.0]},\
               \'yaxis2\': {\'anchor\': \'x2\', \'domain\': [0.0, 1.0]},\
               \'yaxis3\': {\'anchor\': \'x3\', \'domain\': [0.0, 1.0]},\
               \'yaxis4\': {\'anchor\': \'x4\', \'domain\': [0.0, 1.0]}}\
})\
')


class TestModule(unittest.TestCase):
    def test_calculate_centiles(self):
        self.assertEqual(FPTP.calculate_centiles([0,1,2,3,4,5,6,7,8,9,10]),[  9.09,  18.18,  27.27,  36.36,  45.45,  54.55,  63.64,  72.73, 81.82,  90.91, 100.  ])

    def test_check_happy_query_match(self):
        self.assertEqual(FPTP.check_happy_query_match("NA12878", "NA12878"),True)
        self.assertRaises(AssertionError, FPTP.check_happy_query_match, "NA12878" "sample1")

    def test_check_metrics(self):
        self.assertEqual(FPTP.check_metrics("NA12878-query-test.vcf.gz",'AD,DP,GQ,AC,AF,AN,BaseQRankSum'), ['info_AN', 'info_DP', 'info_BaseQRankSum', 'format_GQ', 'format_DP'])
        self.assertEqual(FPTP.check_metrics("NA12878-query-test.vcf.gz",'AD,DP,GQ,AC,AF,AN,BaseQRankSum,MQ,BDA'), ['info_AN', 'info_DP', 'info_BaseQRankSum', 'format_GQ', 'format_DP'])

    def test_check_multiple_query_metrics(self):
        # not implemented yet (multiple samples)        
        pass

    def test_create_plot(self):
        self.assertEqual(str(FPTP.create_plot(test_arrays[0],test_arrays[1])),test_figure)

    def test_decide_bins(self):
        # not implemented yet (multiple samples)        
        pass

    def test_get_args(self):
        pass

    def test_get_output_name(self):
        self.assertEqual(FPTP.get_output_name(['sample1-query.vcf', 'NA12878-query.vcf'], False),"sample1_NA12878_QCDist.html")
        self.assertEqual(FPTP.get_output_name(['NA12878-query-test.vcf.gz']),"TPvsFP_NA12878_QCdist.html")

    def test_get_sample_names(self):
        self.assertEqual(FPTP.get_sample_names('NA12878_test.vcf.gz', 'NA12878-query-test.vcf.gz'), None)

    def test_infer_het_hom(self):
        self.assertEqual(FPTP.infer_het_hom('1/1'), 'homalt')
        self.assertEqual(FPTP.infer_het_hom('0/1'), 'het')
        self.assertEqual(FPTP.infer_het_hom('0/0'), 'homref')

    def test_infer_snp_indel(self):
        self.assertEqual(FPTP.infer_snp_indel('AA', 'A'), 'INDEL')
        self.assertEqual(FPTP.infer_snp_indel('A', 'CAT'), 'INDEL')
        self.assertEqual(FPTP.infer_snp_indel('A', 'G'), 'SNP')

    def test_main(self):
        pass

    def test_make_arrays(self):
        self.assertEqual(FPTP.make_arrays(test_merged_dict, 'info_DP', ['TP','FP'],snp_indel='SNP'), [['TP', 21.0, 79.0, 43.0, 63.0, 64.0, 66.0, 75.0], ['FP']])

    def test_make_html(self):
        pass

    def test_make_plots(self):
        self.assertEqual(str(FPTP.make_plots(test_merged_dict,['info_DP'])),[test_tiled_figure])

    def test_make_report(self):
        # how to test that it writes to file without writing to file?
        pass

    def test_make_tiled_figure(self):
        # need 4 input plots...
        pass

    def test_merge_happy_query(self):
        self.assertEqual(FPTP.merge_happy_query(test_happy_dict, test_query_dict), test_merged_dict)

    def test_merge_samples(self):
        # not implemented yet (multiple samples)        
        pass

    def test_parse_happy(self):
        self.assertEqual(FPTP.parse_happy("NA12878_test.vcf.gz"), test_happy_dict)

    def test_parse_query(self):
        self.assertEqual(FPTP.parse_query("NA12878-query-test.vcf.gz"), test_query_dict)


if __name__ == '__main__':
    unittest.main()