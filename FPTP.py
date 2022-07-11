'''
Script to visualise QC metrics from VCF files as distributions (all variants).
Can compare True Positives to False positives or compare different samples.

Author: Chris Pyatt
'''

# import libraries
import argparse


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
    parser.add_argument('--metrics',
        help='list of metrics to be plotted. The list must be comma-separated with no spaces.'
    )
    args = parser.parse_args()
    # exit gracefully if no arguments given (or missing either file or output)
    if args.file == None or args.output == None:
        parser.print_help()
        sys.exit(1)
    else:
        return args

#parser

#arg happy vcf (optional)

#arg query vcf (at least 1)

#arg list of metrics


def checkMetrics(query, metrics):
    '''
    Takes list of metrics and a query vcf. Checks that all metrics requested are available in query vcf. Returns list of useable metrics for plotting.
    '''
    pass


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

    # generate output
    pass

if __name__ == "__main__":
    main()