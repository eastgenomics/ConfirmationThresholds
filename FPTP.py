'''
Script to visualise QC metrics from VCF files as distributions (all variants).
Can compare True Positives to False positives or compare different samples.

Author: Chris Pyatt
'''

# import libraries
import argparse
import re


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
    # exit gracefully if no arguments given (or missing either file or output)
    if args.query == None or args.output == None:
        parser.print_help()
        sys.exit(1)
    else:
        return args


def checkMetrics(query, metrics):
    '''
    Takes list of metrics and a query VCF. Checks that all metrics requested are available in query VCF. Returns list of useable metrics for plotting.
    '''
    try:
        requestedMetrics = metrics.split(',')
    except:
        print('\nError parsing metrics. Please check formatting.')
        sys.exit(1)
    try:
        with open(query) as file:
            allMetrics = []
            for line in file:
                if line.startswith('##FORMAT') or line.startswith('##INFO'):
                    allMetrics.append(line.split(',')[0].split('=')[-1])
    except:
        print('\nError parsing query VCF. Please check format.')
        sys.exit(1)
    availableMetrics = list(set(requestedMetrics).intersection(allMetrics))
    if args.verbose:
        unavailableMetrics = list(set(requestedMetrics).difference(allMetrics))
        print('\nThe metrics below were requested but not present in the query VCF.\n' + str(unavailableMetrics))
    return availableMetrics


def parseQuery(query):
    '''
    Takes a query vcf filename. Returns a dictionary of relevant metrics and values, paired to variants.
    '''
    pass


def parseHappy(happy):
    '''
    Takes a Hap.py output vcf containing TP & FP calls. Returns a dictionary of variants paired with TP/FP, SNP/INDEL, & het/hom status.
    '''
    pass


def createPlot():
    '''
    Given two arrays of metric values, plot corresponding distributions and return plot object.
    '''
    pass


def getOutputName():
    '''
    Parse input filenames and make output filename
    '''
    pass


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