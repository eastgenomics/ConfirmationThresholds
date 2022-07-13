'''
Script to visualise QC metrics from VCF files as distributions (all variants).
Can compare True Positives to False positives or compare different samples.

Author: Chris Pyatt
'''

# import libraries
import argparse
import re
import sys


def get_args():
    '''
    Parses command line arguments. Returns the arguments as strings.
    '''
    parser = argparse.ArgumentParser(
    
    )
    parser.add_argument(
        '--happy',
        help='filename of hap.py VCF (containing variants and TP/FP calls).'
    )
    parser.add_argument(
        '--verbose',action='store_true',
        help='if enabled, reports which requested metrics are not available in query vcf'
    )
    parser.add_argument(
        '--query',
        help='filename of query VCF (containing variants and metric values). Up to two query sample VCFs may be provided, but if a hap.py VCF is provided via the --happy option, only one (matching) query VCF will be accepted.'
    )
    parser.add_argument('--metrics', default='all',
        help='list of metrics to be plotted. The list must be comma-separated with no spaces. By default, all metrics found in the VCF will be plotted.'
    )
    args = parser.parse_args()
    # exit gracefully if no arguments given or number of input files is not 2
    if (args == None) or (args.query == None) or (len(args.query.split(',')) == 1 and args.happy == None) or (len(args.query.split(',')) == 2 and len(args.happy.split(',')) == 1):
        parser.print_help()
        sys.exit(1)
    else:
        return args


args == None
query == None
query ==  1 and happy == None
query == 2 and happy == 1

def checkHappyQueryMatch():
    happy = args.happy
    query = args.query
    if happy is not None:
        assert happy.split('.')[0] == query.split('.')[0].split('-')[0], f'hap.py and query vcf sample names do not match ({happy} & {query})'
    return


def checkMetrics(query, metrics):
    '''
    Takes list of metrics and a query VCF. Checks that all metrics requested are available in query VCF. Returns list of useable metrics for plotting.
    '''
    try:
        if metrics == 'all':
            requestedMetrics = 'all'
        else:
            requestedMetrics = metrics.split(',')
    except:
        print('\nError parsing metrics. Please check formatting.')
        sys.exit(1)
    try:
        with open(query) as file:
            allInfoMetrics = []
            allFormatMetrics = []
            # parse out metric names - separate info and format metrics in case of identical names (usually DP)
            for line in file:
                if line.startswith('##INFO'):
                    allInfoMetrics.append(line.split(',')[0].split('=')[-1])
                elif line.startswith('##FORMAT'):
                    allFormatMetrics.append(line.split(',')[0].split('=')[-1])
    except:
        print('\nError parsing query VCF. Please check format.')
        sys.exit(1)
    if requestedMetrics == 'all':
        availableInfoMetrics = allInfoMetrics
        availableFormatMetrics = allFormatMetrics
        unavailableMetrics = [[],[]]
    else:
        availableInfoMetrics = list(set(requestedMetrics).intersection(allInfoMetrics))
        availableFormatMetrics = list(set(requestedMetrics).intersection(allFormatMetrics))
        unavailableInfoMetrics = list(set(requestedMetrics).difference(allInfoMetrics))
        unavailableFormatMetrics = list(set(requestedMetrics).difference(allFormatMetrics))
        unavailableMetrics = [unavailableInfoMetrics, unavailableFormatMetrics]
    availableMetrics = [availableInfoMetrics, availableFormatMetrics]
    if args.verbose:
        print('\nThe metrics below were requested but not present in the query VCF.\n' + str(unavailableMetrics))
    return availableMetrics


def parseQuery(query):
    '''
    Takes a query vcf filename. Returns a dictionary of relevant metrics and values, paired to variants.
    '''
    try:
        with open(query) as file:
            variantDict = {}
            for line in file:
                if not line.startswith('#'):
                    variant = line.split('\t')[0] + '_' + line.split('\t')[1] + '_' + line.split('\t')[3] + '_' + line.split('\t')[4]
                    vcf_info = line.split('\t')[7].split(';')
                    vcf_format = line.split('\t')[8].split(':')
                    vcf_genotype = line.split('\t')[9].split(':')
                metricDict = {}
                for item in vcf_info:
                    name = 'info_' + item.split('=')[0]
                    value = item.split('=')[1]
                    metricDict[name] = value
                for i in range(len(vcf_format)):
                    name = 'format_' + vcf_format[i]
                    value = vcf_genotype[i]
                    metricDict[name] = value
                variantDict[variant] = metricDict
    except:
        print('\nError parsing query VCF. Please check format.')
        sys.exit(1)
    return variantDict


def parseHappy(happy):
    '''
    Takes a Hap.py output vcf containing TP & FP calls. Returns a dictionary of variants paired with TP/FP, SNP/INDEL, & het/hom status.
    '''
    try:
        with open(happy) as file:
            variantDict = {}
            for line in file:
                # ignores other info in file e.g. 'CALL_WEIGHT', 'Genotype', 'variant quality for ROC creation', etc.
                if not line.startswith('#'):
                    variant = line.split('\t')[0] + '_' + line.split('\t')[1] + '_' + line.split('\t')[3] + '_' + line.split('\t')[4]
                    queryInfo = line.split('\t')[10].split(':')
                    TPFP = queryInfo[1]
                    SNP_INDEL = queryInfo[5]
                    hethom = queryInfo[6]
                variantDict[variant] = [TPFP,SNP_INDEL,hethom]
    except:
        print('\nError parsing query VCF. Please check format.')
        sys.exit(1)
    return variantDict


def createPlot():
    '''
    Given two arrays of metric values, plot corresponding distributions and return plot object.
    '''
    pass


def getOutputName(file1,file2,happy):
    '''
    Parse input filenames and make output filename.
    '''
    if happy:
        sample1 = file1.split('.')[0]
    else:
        sample1 = file1.split('.')[0].split('-')[0]
    sample2 = file2.split('.')[0].split('-')[0]
    output = f'{sample1}_{sample2}_QCdist.html'
    return output

def makeReport(plots, outFile):
    '''
    Given a list of plot objects, construct an html report and save to file?
    Output name constructed from input filenames.
    '''
    pass


def main():
    # check metrics

    # parse inputs

    # if happy provided, set variable to True, check query length is 1
    # if no happy, set to False, check query length is 2, parse both
    # check metrics for each query & merge lists first

    # generate output
    pass

if __name__ == "__main__":
    main()