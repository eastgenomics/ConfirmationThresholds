'''
Script to visualise QC metrics from VCF files as distributions (all variants).
Can compare True Positives to False positives or compare different samples.

Author: Chris Pyatt
'''

# import libraries
import argparse
import re
import sys


# global variables
VERBOSE = False


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
    # exit gracefully if no arguments given
    if (args == None) or (args.query == None):
        parser.print_help()
        sys.exit(1)
    # set global verbose variable if needed
    if args.verbose:
        VERBOSE = True
    # convert metrics to list
    args.metrics = args.metrics.split(',')
    # convert query to list
    args.query = args.query.split(',')
    # exit gracefully if number of input files is incorrect
    if (len(args.query) == 1 and args.happy == None) or (len(args.query) == 2 and len(args.happy.split(',')) == 1) or (len(args.happy.split(',')) > 1):
        parser.print_help()
        sys.exit(1)
    else:
        return args


def checkHappyQueryMatch(happy, query):
    '''
    Checks that the sample name matches between hap.py VCF and query VCF (as variants will need to be matched between the two to assign TP/FP)
    '''
    if happy is not None:
        assert happy.split('.')[0] == query.split('.')[0].split('-')[0], f'hap.py and query vcf sample names do not match ({happy} & {query})'
    else:
        return True


def checkMultipleQueryMetrics(query, metrics):
    '''
    Takes list of (2) query VCFs and a list of metrics. Returns all shared metrics that are on the list.
    '''
    availableMetrics1 = checkMetrics(query[0], metrics)
    availableMetrics2 = checkMetrics(query[1], metrics)
    availableSharedMetricsInfo = list(set(availableMetrics1[0]).intersection(availableMetrics2[0]))
    availableSharedMetricsFormat = list(set(availableMetrics1[1]).intersection(availableMetrics2[1]))
    availableSharedMetrics = [availableSharedMetricsInfo, availableSharedMetricsFormat]
    return availableSharedMetrics


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
    if VERBOSE:
        print(f'\nThe metrics below were requested but not present in the query VCF ({query}).\n' + str(unavailableMetrics))
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


def createPlot(array1, array2):
    '''
    Given two arrays of metric values, plot corresponding distributions and return plot object.
    '''
    pass


def getOutputName(files, happy=True):
    '''
    Parse input filenames and make output filename. Takes list of query files (1 or more) and boolean to toggle whether matching happy files are inputs. Returns sample1_sample2_QCdist.html if two sample inputs, otherwise TPvsFP_samplename_QCdist.html for happy inputs.
    '''
    if happy:
        sample1 = 'TPvsFP'
    else:
        sample1 = files[0].split('.')[0].split('-')[0]
    sample2 = files[1].split('.')[0].split('-')[0]
    output = f'{sample1}_{sample2}_QCdist.html'
    return output


def makeReport(plots, outFile):
    '''
    Given a list of plot objects, construct an html report and save to file?
    Output name constructed from input filenames.
    '''
    pass


def mergeSamples(sample1, sample2):
    '''
    Take two dictionaries, merge them, return merged dictionary.
    '''
    mergedDict = sample1 | sample2
    return mergedDict


def makeArrays():
    '''
    Take 
    '''
    pass


def main():
    # get command line arguments
    args = get_args()
    if args.happy:
        # check that sample names match between happy and query files
        checkHappyQueryMatch(args.happy, args.query[0])
        # check metrics are available
        metrics = checkMetrics(args.query[0], args.metrics)
        # parse inputs
        sample1 = parseHappy(args.happy)
        sample2 = parseQuery(args.query[0])
    else:
        # check metrics are available
        metrics = checkMultipleQueryMetrics(args.query, args.metrics)
        # parse inputs
        sample1 = parseQuery(args.query[0])
        sample2 = parseQuery(args.query[1])

    # list of plot objects to insert into report
    plots = []
    for metric in metrics:
        # list of 4 plots
        metric_plots = []
        # plot snps & append to metric_plots
        # plot indels & append to metric_plots
        # plot hets & append to metric_plots
        # plot homs & append to metric_plots
        plots.append(metric_plots)
    
    # for each metric in plots list, add to report ????????

    # generate output
    if args.happy:
        output = getOutputName(args.query)
    else:
        output = getOutputName(args.query, False)


if __name__ == "__main__":
    main()