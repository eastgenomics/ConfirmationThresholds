'''
Script to visualise QC metrics from VCF files as distributions (all variants).
Can compare True Positives to False positives or compare different samples.

Author: Chris Pyatt
'''

# import libraries
import argparse
from msilib.schema import Error
from operator import mod
import re
import sys
import pandas as pd
import numpy as np
#import plotly.plotly as py
import plotly.figure_factory as ff
import plotly.graph_objects as go
import gzip
from numpy.linalg import inv, det
from IPython.display import HTML


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
                    variant = (
                        line.split('\t')[0] + '_' + line.split('\t')[1] + '_'
                        + line.split('\t')[3] + '_' + line.split('\t')[4]
                    )
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
                    for i in range(len(vcf_format)):
                        name = 'format_' + vcf_format[i]
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
                    variant = (
                        line.split('\t')[0] + '_' + line.split('\t')[1] + '_'
                        + line.split('\t')[3] + '_' + line.split('\t')[4]
                    )
                    query_info = line.split('\t')[10].split(':')
                    cat_dict['TPFP_or_samplename'] = query_info[1]
                    cat_dict['snp_indel'] = query_info[5]
                    cat_dict['HETHOM'] = query_info[6].rstrip()
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


def create_plot(array1, array2, name):
    '''
    Given two arrays of metric values, plot corresponding distributions and
    return plot object. Name variable to be title of plot?
    '''
    labels = [array1.pop(0), array2.pop(0)]
    if len(array1) < 1 or len(array2) < 1:
        # do something to indicate insufficient data for this metric combo??
        return None
    # convert arrays to dataframe with column headers
    plot_df = pd.DataFrame(
        {
            labels[0]: np.random.randn(200),
            labels[1]: np.random.randn(200)+1
        }
        )
    # make distribution plot object
    fig = ff.create_distplot(
        [plot_df[c] for c in plot_df.columns], plot_df.columns, bin_size=.25
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


def make_report(plots, outFile):
    '''
    Given a list of plot objects, construct an html report and save to file?
    Output name constructed from input filenames.
    '''
    pass
    #for plot_obj in plots:
    #    plot_url = py.plot(plot_obj, filename='plot_object', auto_open=False,)
    #    print(plot_url)
    #return


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
        fig = make_tiled_figure([snp_plot, indel_plot, het_plot, hom_plot])
        plot_list.append(fig)
    return plot_list


def make_tiled_figure(subfigs):
    '''
    Take list of figures ( figure factory plot objects) to be combined into
    tiled image. Return single figure object with tiled subplots.
    '''
    for i in range(len(subfigs)):
        # check subfigure is not empty
        if not subfigs[i]:
            continue
        # subplot suffixes cannot be 0 indexed
        index = i + 1
        # initialize xaxis and yaxis
        subfigs[i]['layout'][f'xaxis{index}'] = {}
        subfigs[i]['layout'][f'yaxis{index}'] = {}

        for j in range(len(subfigs[i].data)):
            subfigs[i].data[j].xaxis=f'x{index}'
            subfigs[i].data[j].yaxis=f'y{index}'

        subfigs[i]['layout'][f'xaxis{index}'].update({'anchor': f'y{index}'})
        subfigs[i]['layout'][f'yaxis{index}'].update({'anchor': f'x{index}', 'domain': [(.25*index), 1-(.25*index)]})

    fig = go.Figure()

    for subfig in subfigs:
        if subfig:
            fig.add_traces(subfig.data[0])
            fig.layout.update(subfig.layout)
    return fig


def make_html(fig):
    '''
    Temporary function to help me visualise the output before making a proper html file.
    '''
    import dash
    from dash import dcc
    from dash import html

    app = dash.Dash()
    app.layout = html.Div([
        dcc.Graph(figure=fig)
    ])

    app.run_server(debug=True, use_reloader=False)  # Turn off reloader if inside Jupyter


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

    # for each metric in plots list, add to report ????????
    make_report(plots, output)

    print(f'Type of plots: {type(plots)}')
    make_html(plots[0])

    #for i in plots:
    #    print(f'Type of individual: {type(i)}')
    #    make_html(i)

if __name__ == "__main__":
    main()