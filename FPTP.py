'''
Script to visualise QC metrics from VCF files as distributions (all variants).
Can compare True Positives to False positives or compare different samples.

Author: Chris Pyatt
'''

# import libraries
import argparse
import re
import sys
import pandas as pd
import plotly
import gzip


# global variables
VERBOSE = False
# names used several times so let's just save them here
SAMPLE1_NAME = ''
SAMPLE2_NAME = ''


def get_args():
    '''
    Parses command line arguments. Returns the arguments as strings.
    '''
    parser = argparse.ArgumentParser(
    
    )
    parser.add_argument(
        '--happy',
        help='filename of hap.py VCF (containing variants and TP/FP calls). Gzipped.'
    )
    parser.add_argument(
        '--verbose',action='store_true',
        help='if enabled, reports which requested metrics are not available in query vcf'
    )
    parser.add_argument(
        '--query',
        help='filename of query VCF (containing variants and metric values). Gzipped. Up to two query sample VCFs may be provided, but if a hap.py VCF is provided via the --happy option, only one (matching) query VCF will be accepted. Input data must be normalised VCF.'
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
        global VERBOSE
        VERBOSE = True
    # convert metrics to list
    #args.metrics = args.metrics.split(',')
    # convert query to list
    args.query = args.query.split(',')
    # exit gracefully if number of input files is incorrect
    if (len(args.query) == 1 and args.happy == None) or (len(args.query) == 2 and len(args.happy.split(',')) == 1) or (len(args.happy.split(',')) > 1):
        parser.print_help()
        sys.exit(1)
    else:
        return args


def getSampleNames(sample1, sample2):
    try:
        # is declaring as global necessary if declared at top of script?
        global SAMPLE1_NAME
        global SAMPLE2_NAME
        SAMPLE1_NAME = re.split('[.\-]', sample1)[0]
        SAMPLE2_NAME = re.split('[.\-]', sample2)[0]
    except:
        print(f'One of {sample1} or {sample2} does not match the expected pattern. Expected something like samplename-metadata.vcf.gz')
        sys.exit(1)


def checkHappyQueryMatch(happy, query):
    '''
    Checks that the sample name matches between hap.py VCF and query VCF (as variants will need to be matched between the two to assign TP/FP)
    '''
    if happy:
        assert happy == query, f'hap.py and query vcf sample names do not match ({happy} & {query})'
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
        with gzip.open(query, 'rt') as file:
            allInfoMetrics = []
            allFormatMetrics = []
            # parse out metric names - separate info and format metrics in case of identical names (usually DP)
            for line in file:
                if line.startswith('##INFO'):
                    allInfoMetrics.append(line.split(',')[0].split('=')[-1])
                elif line.startswith('##FORMAT'):
                    allFormatMetrics.append(line.split(',')[0].split('=')[-1]) # add stop after header to prevent parsing whole file here
    except:
        print('\nError retreiving query VCF metrics. Please check format.')
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
        with gzip.open(query, 'rt') as file:
            variantDict = {}
            for line in file:
                if not line.startswith('#'):
                    variant = line.split('\t')[0] + '_' + line.split('\t')[1] + '_' + line.split('\t')[3] + '_' + line.split('\t')[4]
                    vcf_info = line.split('\t')[7].split(';')
                    vcf_format = line.split('\t')[8].split(':')
                    vcf_genotype = line.split('\t')[9].split(':')
                    metricDict = {}
                    # get samplename
                    metricDict['TPFP_or_samplename'] = query.split('.')[0].split('-')[0]
                    # infer SNP/INDEL status based on ref/alt fields
                    metricDict['SNP_INDEL'] = inferSnpIndel(variant.split('_')[2], variant.split('_')[3])
                    # infer het/hom status based on genotype field
                    metricDict['HETHOM'] = inferHetHom(vcf_genotype[0])
                    # add info metrics to dictionary
                    for item in vcf_info:
                        name = 'info_' + item.split('=')[0]
                        # catch flag metrics
                        try:
                            value = item.split('=')[1]
                        except IndexError:
                            value = True
                        metricDict[name] = value
                    # add format/genotype metrics to dictionary
                    for i in range(len(vcf_format)):
                        name = 'format_' + vcf_format[i]
                        value = vcf_genotype[i]
                        metricDict[name] = value
                    variantDict[variant] = metricDict
    except:
        print('\nError parsing query VCF. Please check format.')
        sys.exit(1)
    return variantDict


def inferHetHom(genotype):
    '''
    Take genotype field, return het or homalt label
    '''
    if genotype == '0/0':
        hethom = 'homref'
    elif genotype == '0/1' or genotype == '1/0':
        hethom = 'het'
    elif genotype == '1/1':
        hethom = 'homalt'
    else:
        # not sure whether exiting at this point is best given could be one of thousands - maybe save as N/A and deal with it elsewhere?
        print(f'Genotype {genotype} not recognised. Het/hom cannot be inferred.')
        sys.exit(1)
    return hethom


def inferSnpIndel(ref, alt):
    '''
    Take ref & alt fields, return SNP or INDEL label. Input VCF must be normalised for this method to work (i.e. one variant per position).
    '''
    if ',' in alt:
        print(f'Comma present in alt field ({alt}), suggesting non-normalised VCF. Please provide normalised VCF data.')
        sys.exit(1)
    elif len(ref) == 1 or len(alt) == 1:
        SNP_INDEL = 'SNP'
    else:
        SNP_INDEL = 'INDEL'
    return SNP_INDEL


def parseHappy(happy):
    '''
    Takes a Hap.py output vcf containing TP & FP calls. Returns a dictionary of variants paired with TP/FP, SNP/INDEL, & het/hom status.
    '''
    try:
        with gzip.open(happy, 'rt') as file:
            variantDict = {}
            for line in file:
                # ignores other info in file e.g. 'CALL_WEIGHT', 'Genotype', 'variant quality for ROC creation', etc.
                catDict = {}
                if not line.startswith('#'):
                    variant = line.split('\t')[0] + '_' + line.split('\t')[1] + '_' + line.split('\t')[3] + '_' + line.split('\t')[4]
                    queryInfo = line.split('\t')[10].split(':')
                    catDict['TPFP_or_samplename'] = queryInfo[1]
                    catDict['SNP_INDEL'] = queryInfo[5]
                    catDict['HETHOM'] = queryInfo[6].rstrip()
                else:
                    continue
                variantDict[variant] = catDict
    except:
        print('\nError parsing happy VCF. Please check format.')
        sys.exit(1)
    return variantDict


def createPlot(array1, array2):
    '''
    Given two arrays of metric values, plot corresponding distributions and return plot object.
    '''
    labels = [array1.pop(0), array2.pop(0)]
    fig = ff.create_distplot(variables, labels, show_hist=False)
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


def mergeHappyQuery(happy, query):
    '''
    Take happy dict (output of parseHappy) and query dict (output of parseQuery) & return new dictionary with matched variants and merged values
    '''
    mergedDict = {}
    missingVariants = []
    for variant in happy:
        try:
            mergedDict[variant] = {**happy[variant], **query[variant]}
        except KeyError:
            missingVariants.append(variant)
    if VERBOSE:
        print(f'\nThe variants below were present in the hap.py VCF but not the query VCF.\n' + str(missingVariants))
    return mergedDict


def mergeSamples(sample1, sample2):
    '''
    Take two dictionaries, merge them, return merged dictionary. Keys now contain samplename in case samples share variants.
    '''
    mergedDict = {}
    for variant in sample1:
        new_key = SAMPLE1_NAME + '-' + variant
        mergedDict[new_key] = sample1[variant]
    for variant in sample2:
        new_key = SAMPLE2_NAME + '-' + variant
        mergedDict[new_key] = sample2[variant]
    return mergedDict


def makeArrays(data, metric, fptp, snp_indel=None, hethom=None):
    '''
    Take merged dictionary, return two arrays (happy vs query, or sample1 vs sample2) of relevant metric values for plotting, split according to category (SNP/INDEL, het/hom).
    '''
    array1 = [fptp[0]]
    array2 = [fptp[1]]
    if snp_indel:
        filtered_keys_1 = [k for k,v in data.items() if v['TPFP_or_samplename'] == fptp[0] and v['SNP_INDEL'] == snp_indel]
        filtered_keys_2 = [k for k,v in data.items() if v['TPFP_or_samplename'] == fptp[1] and v['SNP_INDEL'] == snp_indel]
    if hethom:
        filtered_keys_1 = [k for k,v in data.items() if v['TPFP_or_samplename'] == fptp[0] and v['HETHOM'] == hethom]
        filtered_keys_2 = [k for k,v in data.items() if v['TPFP_or_samplename'] == fptp[1] and v['HETHOM'] == hethom]
    else:
        filtered_keys_1 = [k for k,v in data.items() if v['TPFP_or_samplename'] == fptp[0]]
        filtered_keys_2 = [k for k,v in data.items() if v['TPFP_or_samplename'] == fptp[1]]
    for item in filtered_keys_1:
        array1.append(data[item][metric])
    for item in filtered_keys_2:
        array2.append(item[metric])
    return [array1, array2]


def makeArrayList(data, metrics, happy=True):
    if happy:
        fptp = ['TP', 'FP']
    else:
        fptp = [SAMPLE1_NAME, SAMPLE2_NAME]
    for metric in metrics:
        snp_arrays = makeArrays(data, metric, fptp, snp_indel='SNP')
        indel_arrays = makeArrays(data, metric, fptp, snp_indel='INDEL')
        het_arrays = makeArrays(data, metric, fptp, hethom='het')
        hom_arrays = makeArrays(data, metric, fptp, hethom='homalt')
    return [snp_arrays, indel_arrays, het_arrays, hom_arrays]


def main():
    # get command line arguments
    args = get_args()
    if args.happy:
        getSampleNames(args.happy, args.query[0])
        # check that sample names match between happy and query files
        checkHappyQueryMatch(SAMPLE1_NAME, SAMPLE2_NAME)
        # check metrics are available
        metrics = checkMetrics(args.query[0], args.metrics)
        # parse inputs
        sample1 = parseHappy(args.happy)
        sample2 = parseQuery(args.query[0])
        # merge input dicts
        merged_data = mergeHappyQuery(sample1, sample2)
        # make 4 arrays for snp, indel, het, hom plots
        arrays = makeArrayList(merged_data, metrics)
    else:
        getSampleNames(args.query[0], args.query[1])
        # check metrics are available
        print("Dual sample inputs not yet implemented\n")
        sys.exit(1)
        metrics = checkMultipleQueryMetrics(args.query, args.metrics)
        # parse inputs
        sample1 = parseQuery(args.query[0])
        sample2 = parseQuery(args.query[1])
        # merge input dicts
        merged_data = mergeSamples(sample1, sample2, False)
        # make 4 arrays for snp, indel, het, hom plots
        arrays = makeArrayList(merged_data, metrics, happy=False)


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
        output = getOutputName([args.happy, args.query[0]])
    else:
        output = getOutputName(args.query, False)


if __name__ == "__main__":
    main()