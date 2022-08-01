'''
Script to visualise QC metrics from VCF files as distributions (all variants).
Can compare True Positives to False positives or compare different samples.

Author: Chris Pyatt
'''

# import libraries
import argparse
import re
import sys
import math
import numpy as np
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
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
        help=(
            'filename of hap.py VCF (containing variants and TP/FP calls). \
                Gzipped.'
        )
    )
    parser.add_argument(
        '--verbose', action='store_true',
        help=(
            'if enabled, reports which requested metrics are not available \
                in query vcf'
        )
    )
    parser.add_argument(
        '--query', required=True,
        help=(
            'filename of query VCF (containing variants and metric values). \
                Gzipped. Up to two query sample VCFs may be provided, but if \
                    a hap.py VCF is provided via the --happy option, only one \
                        (matching) query VCF will be accepted. Input data \
                            must be normalised VCF.'
        )
    )
    parser.add_argument(
        '--metrics', default='all',
        help=(
            'list of metrics to be plotted. The list must be comma-separated \
                with no spaces. By default, all metrics found in the VCF will \
                    be plotted.'
        )
    )
    args = parser.parse_args()
    # exit gracefully if no arguments given
    if not args or not args.query:
        parser.print_help()
        sys.exit(1)
    # set global verbose variable if needed
    if args.verbose:
        global VERBOSE
        VERBOSE = True
    # convert query to list
    args.query = args.query.split(',')
    # exit gracefully if number of input files is incorrect
    if (
        (len(args.query) == 1 and not args.happy)
        or (len(args.query) == 2 and len(args.happy.split(',')) == 1)
        or (len(args.happy.split(',')) > 1)
    ):
        parser.print_help()
        sys.exit(1)
    else:
        return args


def get_sample_names(sample1, sample2):
    '''
    Takes two input files (should be in the format \
    samplename-some-metadata-fields.vcf.gz) and returns just \
    the sample names as global variables
    '''
    try:
        # is declaring as global necessary if declared at top of script?
        global SAMPLE1_NAME
        global SAMPLE2_NAME
        SAMPLE1_NAME = re.split('[.-]', sample1)[0]
        SAMPLE2_NAME = re.split('[.-]', sample2)[0]
    except Exception as error:
        print(
            f'One of {sample1} or {sample2} does not match the expected \
                pattern. Expected something like samplename-metadata.vcf.gz'
            f'Error: {error}'
            )
        sys.exit(1)


def check_happy_query_match(happy, query):
    '''
    Checks that the sample name matches between hap.py VCF and query VCF
    (as variants will need to be matched between the two to assign TP/FP)
    '''
    if happy:
        assert happy == query, f'hap.py and query vcf sample names do not \
            match ({happy} & {query})'
    else:
        return True


def check_multiple_query_metrics(query, metrics):
    '''
    Takes list of (2) query VCFs and a list of metrics. Returns all shared
    metrics that are on the list.
    '''
    available_metrics1 = check_metrics(query[0], metrics)
    available_metrics2 = check_metrics(query[1], metrics)
    available_shared_metrics_info = list(
        set(available_metrics1[0]).intersection(available_metrics2[0])
        )
    available_shared_metrics_format = list(
        set(available_metrics1[1]).intersection(available_metrics2[1])
        )
    available_shared_metrics = [available_shared_metrics_info,
                                available_shared_metrics_format]
    return available_shared_metrics


def check_metrics(query, metrics):
    '''
    Takes list of metrics and a query VCF. Checks that all metrics requested
    are available in query VCF. Returns list of useable metrics for plotting.
    '''
    try:
        if metrics == 'all':
            requested_metrics = 'all'
        else:
            requested_metrics = metrics.split(',')
    except Exception as error:
        print(
            '\nError parsing metrics. Please check formatting.'
            f'Error: {error}'
            )
        sys.exit(1)
    try:
        with gzip.open(query, 'rt') as file:
            all_info_metrics = []
            all_format_metrics = []
            # parse out metric names - separate info and format metrics in case of identical names (usually DP)
            for line in file:
                # grab only metrics where the type is float or integer (as others cannot be plotted) and also where the number of values is constrained to 1. The latter constraint will be modified to accept lists of values in later versions of this tool.
                if line.startswith('##INFO') and ('Integer' in line.split(',')[2] or 'Float' in line.split(',')[2]) and line.split(',')[1] == 'Number=1':
                    all_info_metrics.append(line.split(',')[0].split('=')[-1])
                elif line.startswith('##FORMAT') and ('Integer' in line.split(',')[2] or 'Float' in line.split(',')[2]) and line.split(',')[1] == 'Number=1':
                    all_format_metrics.append(line.split(',')[0].split('=')[-1])
                elif line.startswith('#CHROM'):
                    break
    except Exception as error:
        print(
            '\nError retreiving query VCF metrics. Please check format.'
            f'Error: {error}'
            )
        sys.exit(1)
    if requested_metrics == 'all':
        available_info_metrics = all_info_metrics
        available_format_metrics = all_format_metrics
        unavailable_metrics = [[], []]
    else:
        available_info_metrics = list(
            set(requested_metrics).intersection(all_info_metrics)
            )
        available_format_metrics = list(
            set(requested_metrics).intersection(all_format_metrics)
            )
        unavailable_info_metrics = list(
            set(requested_metrics).difference(all_info_metrics)
            )
        unavailable_format_metrics = list(
            set(requested_metrics).difference(all_format_metrics)
            )
        unavailable_metrics = [unavailable_info_metrics,
                               unavailable_format_metrics]
    # add prefixes to metrics list to distinguish duplicate metric names in \
    # INFO and FORMAT fields
    available_metrics = []
    for metric in available_info_metrics:
        available_metrics.append(f'info_{metric}')
    for metric in available_format_metrics:
        available_metrics.append(f'format_{metric}')
    if VERBOSE:
        print(f'\nThe metrics below were requested but not present in the \
            query VCF ({query}).\n' + str(unavailable_metrics))
    return available_metrics


def parse_query(query, happy=True):
    '''
    Takes a query vcf filename. Returns a dictionary of relevant metrics and
    values, paired to variants.
    '''
    try:
        with gzip.open(query, 'rt') as file:
            variant_dict = {}
            for line in file:
                if not line.startswith('#'):
                    chrom = line.split('\t')[0]
                    pos = line.split('\t')[1]
                    ref = line.split('\t')[3]
                    alt = line.split('\t')[4]
                    variant = (f'{chrom}_{pos}_{ref}_{alt}')
                    vcf_info = line.split('\t')[7].split(';')
                    vcf_format = line.split('\t')[8].split(':')
                    vcf_genotype = line.rstrip().split('\t')[9].split(':')
                    metric_dict = {}
                    if not happy:
                        # get samplename
                        metric_dict['TPFP_or_samplename'] = (
                                        query.split('.')[0].split('-')[0]
                        )
                        # infer SNP/INDEL status based on ref/alt fields
                        metric_dict['snp_indel'] = infer_snp_indel(
                            variant.split('_')[2], variant.split('_')[3]
                            )
                        # infer het/hom status based on genotype field
                        metric_dict['HETHOM'] = infer_het_hom(vcf_genotype[0])
                    # add info metrics to dictionary
                    for item in vcf_info:
                        name = 'info_' + item.split('=')[0]
                        # catch flag metrics
                        try:
                            value = item.split('=')[1]
                        except IndexError:
                            value = True
                        metric_dict[name] = value
                    # add format/genotype metrics to dictionary
                    for i, metric in enumerate(vcf_format):
                        name = 'format_' + metric
                        value = vcf_genotype[i]
                        metric_dict[name] = value
                    variant_dict[variant] = metric_dict
    except Exception as error:
        print(
            '\nError parsing query VCF. Please check format.'
            f'Error: {error}'
            )
        sys.exit(1)
    return variant_dict


def infer_het_hom(genotype):
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
        # not sure whether exiting at this point is best given could be one of\
        # thousands - maybe save as N/A and deal with it elsewhere?
        print(
            f'Genotype {genotype} not recognised. Het/hom cannot be inferred.'
            )
        sys.exit(1)
    return hethom


def infer_snp_indel(ref, alt):
    '''
    Take ref & alt fields, return SNP or INDEL label. Input VCF must be
    normalised for this method to work (i.e. one variant per position).
    '''
    if ',' in alt:
        print(
            f'Comma present in alt field ({alt}), suggesting non-normalised \
            VCF. Please provide normalised VCF data.'
            )
        sys.exit(1)
    elif len(ref) == 1 or len(alt) == 1:
        snp_indel = 'SNP'
    else:
        snp_indel = 'INDEL'
    return snp_indel


def parse_happy(happy):
    '''
    Takes a Hap.py output vcf containing TP & FP calls. Returns a dictionary
    of variants paired with TP/FP, SNP/INDEL, & het/hom status.
    '''
    try:
        with gzip.open(happy, 'rt') as file:
            variant_dict = {}
            for line in file:
                # ignores other info in file e.g. 'CALL_WEIGHT', 'Genotype', \
                # 'variant quality for ROC creation', etc.
                cat_dict = {}
                if not line.startswith('#'):
                    chrom = line.split('\t')[0]
                    pos = line.split('\t')[1]
                    ref = line.split('\t')[3]
                    alt = line.split('\t')[4]
                    variant = (f'{chrom}_{pos}_{ref}_{alt}')
                    query_format = line.split('\t')[8].split(':')
                    query_values = line.rstrip().split('\t')[10].split(':')
                    # add format/query metrics to dictionary (format varies)
                    # so can't rely on simple string splitting
                    query_dict = {}
                    for i, metric in enumerate(query_format):
                        value = query_values[i]
                        query_dict[metric] = value
                    cat_dict['TPFP_or_samplename'] = query_dict['BD']
                    cat_dict['snp_indel'] = query_dict['BVT']
                    cat_dict['HETHOM'] = query_dict['BLT']
                else:
                    continue
                variant_dict[variant] = cat_dict
    except Exception as error:
        print(
            '\nError parsing happy VCF. Please check format.'
            f'Error: {error}'
            )
        sys.exit(1)
    return variant_dict


def decide_bins(array):
    '''
    Take numpy array & use number and range of values to determine
    appropriate bin size for histopgram. Returns bin size as integer.
    '''
    array_length = len(array)
    array_range = max(array) - min(array)
    # Use modified Sturge's rule to decide number of bins
    num_bins = (1 + 3.322 * math.log10(array_length)) * 3
    # Divide the range of values by the number of bins to get bin size
    bin_size = array_range / num_bins
    return bin_size


def create_plot(array1, array2, name):
    '''
    Given two arrays of metric values, plot corresponding distributions and
    return plot object. Name variable to be title of plot?
    '''
    labels = [array1.pop(0), array2.pop(0)]
    colours = ['#eb8909', '#09ebeb']
    if len(array1) < 1 or len(array2) < 1:
        # do something to indicate insufficient data for this metric combo??
        return None
    # decide bin sizes based on array1 (should be TPs so the longer dataset)
    bin_size = decide_bins(array1)
    # convert arrays to dataframe with column headers
    hist_data = [np.array(array1), np.array(array2)]
    # make distribution plot object - no curves as gets broken by symmetrical
    # matrix (all values the same in this case)
    fig = ff.create_distplot(
                             hist_data, labels, bin_size=bin_size,
                             colors=colours, show_curve=False
                             )
    return fig


def get_output_name(files, happy=True):
    '''
    Parse input filenames and make output filename. Takes list of query files
    (1 or more) and boolean to toggle whether matching happy files are inputs.
    Returns sample1_sample2_QCdist.html if two sample inputs, otherwise
    TPvsFP_samplename_QCdist.html for happy inputs.
    '''
    if happy:
        sample1 = 'TPvsFP'
    else:
        sample1 = files[0].split('.')[0].split('-')[0]
    sample2 = files[1].split('.')[0].split('-')[0]
    output = f'{sample1}_{sample2}_QCdist.html'
    return output


def make_html(plots):
    '''
    Given a list of plot objects, construct an html report and save to file?
    Output name constructed from input filenames.
    '''
    css = 'https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css'
    html_string = (
                    f'<html><head><link rel="stylesheet" href="{css}">\
                    <style>body{{ margin:0 100; background:whitesmoke; }}\
                    </style></head><body>'
                )
    for i in plots:
        html_string = html_string + '<h2>Section 1: Apple Inc. (AAPL) stock in 2014</h2>'
        html_string = html_string + (
                        f'<iframe width="1000" height="550" frameborder="0" \
                        seamless="seamless" scrolling="no" \
                        src="{i}.embed?width=800&height=550"></iframe>'
                       )
    html_string = html_string + '</body></html>'
    return html_string


def make_report(html_string, outfile):
    '''
    Take html string (output of make_html) and save to output file
    '''
    try:
        with open(outfile, 'w') as out:
            out.write(html_string)
        return 0
    except Exception as error:
        print(
            '\nError saving output to html file.'
            f'Error: {error}'
            )
        sys.exit(1)


def merge_happy_query(happy, query):
    '''
    Take happy dict (output of parse_happy) and query dict (output of
    parse_query) & return new dictionary with matched variants and merged
    values
    '''
    merged_dict = {}
    missing_variants = []
    for variant in happy:
        try:
            merged_dict[variant] = {**happy[variant], **query[variant]}
        except KeyError:
            missing_variants.append(variant)
    if VERBOSE:
        print(
            f'\nThe variants below were present in the hap.py VCF but not the \
                query VCF:\n {str(missing_variants)}'
            )
    return merged_dict


def merge_samples(sample1, sample2):
    '''
    Take two dictionaries, merge them, return merged dictionary. Keys now
    contain samplename in case samples share variants.
    '''
    merged_dict = {}
    for variant in sample1:
        new_key = SAMPLE1_NAME + '-' + variant
        merged_dict[new_key] = sample1[variant]
    for variant in sample2:
        new_key = SAMPLE2_NAME + '-' + variant
        merged_dict[new_key] = sample2[variant]
    return merged_dict


def make_arrays(data, metric, fptp, snp_indel=None, hethom=None):
    '''
    Take merged dictionary, return two arrays (happy vs query, or sample1 vs
    sample2) of relevant metric values for plotting, split according to
    category (SNP/INDEL, het/hom).
    '''
    array1 = [fptp[0]]
    array2 = [fptp[1]]
    if snp_indel:
        filtered_keys_1 = (
            [k for k, v in data.items() if v['TPFP_or_samplename'] == fptp[0]
             and v['snp_indel'] == snp_indel]
        )
        filtered_keys_2 = (
            [k for k, v in data.items() if v['TPFP_or_samplename'] == fptp[1]
             and v['snp_indel'] == snp_indel]
        )
    if hethom:
        filtered_keys_1 = (
            [k for k, v in data.items() if v['TPFP_or_samplename'] == fptp[0]
             and v['HETHOM'] == hethom]
        )
        filtered_keys_2 = (
            [k for k, v in data.items() if v['TPFP_or_samplename'] == fptp[1]
             and v['HETHOM'] == hethom]
        )
    else:
        filtered_keys_1 = (
            [k for k, v in data.items() if v['TPFP_or_samplename'] == fptp[0]]
        )
        filtered_keys_2 = (
            [k for k, v in data.items() if v['TPFP_or_samplename'] == fptp[1]]
        )
    for item in filtered_keys_1:
        # try to append metric value to array but catch occurrences where \
        # metric is not present for that variant (usually metrics like \
        # BaseQRankSum, ClippingRankSum, ExcessHet, etc.)
        try:
            array1.append(float(data[item][metric]))
        except KeyError:
            pass
    for item in filtered_keys_2:
        try:
            array2.append(float(data[item][metric]))
        except KeyError:
            pass
    return [array1, array2]


def make_plots(data, metrics, happy=True):
    '''
    Take merged data and list of metrics. Return list of plot objects.
    '''
    plot_list = []
    if happy:
        fptp = ['TP', 'FP']
    else:
        fptp = [SAMPLE1_NAME, SAMPLE2_NAME]
    for metric in metrics:
        # make filtered arrays for each variant category
        snp_arrays = make_arrays(data, metric, fptp, snp_indel='SNP')
        indel_arrays = make_arrays(data, metric, fptp, snp_indel='INDEL')
        het_arrays = make_arrays(data, metric, fptp, hethom='het')
        hom_arrays = make_arrays(data, metric, fptp, hethom='homalt')
        # make plots for each category
        snp_plot = create_plot(snp_arrays[0], snp_arrays[1], 'SNP')
        indel_plot = create_plot(indel_arrays[0], indel_arrays[1], 'INDEL')
        het_plot = create_plot(het_arrays[0], het_arrays[1], 'HET')
        hom_plot = create_plot(hom_arrays[0], hom_arrays[1], 'HOM')
        # make tiled figure with all of the above
        fig = make_tiled_figure(
                                [snp_plot, indel_plot, het_plot, hom_plot],
                                metric
                                )
        plot_list.append(fig)
    return plot_list


def make_tiled_figure(subfigs, metric):
    '''
    Take list of figures ( figure factory plot objects) to be combined into
    tiled image. Return single figure object with tiled subplots.
    '''
    fig = make_subplots(rows=2, cols=2)
    # decide on position and add subfigures to plot
    for i, subfig in enumerate(subfigs):
        if i in (1, 2):
            row_val = 1
        else:
            row_val = 2
        if i in (1, 3):
            col_val = 1
        else:
            col_val = 2
        if subfig:
            fig.add_trace(subfig.data[0], row=row_val, col=col_val)
    # specify plot size and title
    fig.update_layout(height=1000, width=1000, title_text=metric)
    return fig


def main():
    '''
    Main app code calling other functions.
    '''
    # get command line arguments
    args = get_args()
    if args.happy:
        get_sample_names(args.happy, args.query[0])
        # check that sample names match between happy and query files
        check_happy_query_match(SAMPLE1_NAME, SAMPLE2_NAME)
        # check metrics are available
        metrics = check_metrics(args.query[0], args.metrics)
        # parse inputs
        sample1 = parse_happy(args.happy)
        sample2 = parse_query(args.query[0])
        # merge input dicts
        merged_data = merge_happy_query(sample1, sample2)
        # make 4 arrays for snp, indel, het, hom plots, for each metric & make plot objects
        plots = make_plots(merged_data, metrics)
    else:
        get_sample_names(args.query[0], args.query[1])
        # check metrics are available
        print("Dual sample inputs not yet implemented\n")
        sys.exit(1)
        metrics = check_multiple_query_metrics(args.query, args.metrics)
        # parse inputs
        sample1 = parse_query(args.query[0], False)
        sample2 = parse_query(args.query[1], False)
        # merge input dicts
        merged_data = merge_samples(sample1, sample2, False)
        # make 4 arrays for snp, indel, het, hom plots
        plots = make_plots(merged_data, metrics, happy=False)
    
    # generate output filename
    if args.happy:
        output = get_output_name([args.happy, args.query[0]])
    else:
        output = get_output_name(args.query, False)

    # for each metric in plots list, add to report
    html_output = make_html(plots)
    make_report(html_output, output)


if __name__ == "__main__":
    main()
