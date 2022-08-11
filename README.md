# ConfirmationThresholds
A tool to plot distributions of VCF metrics separated by True Positive / False Positive

Use the output of this tool to decide QC thresholds above or below which Sanger confirmation may not be necessary.

USAGE:
python3 FPTP.py --happy [happy vcf] --query [query vcf] [options]


INPUT:
1) a hap.py VCF (found in the zipped folder created by vcfeval_hap.py).
2) a matching query VCF (the Sentieon output used to generate the hap.py VCF).
*these files must match in order for the tool's output to be informative*

OPTIONAL INPUTS:
1) a list of metrics to plot. By default all metrics will be plotted but you can choose a subset if desired. Any requested metric that is not found in the input VCF will be ignored.
2) a verbose flag. If enabled, the tool will print extra information to the terminal, including the requested metrics that were not available in the input VCF.

OUTPUT:
1) an HTML report with interactive plots for all metrics. Each metric will, if possible, have four plots with data split according to SNP/INDEL and HET/HOM status. The text revealed by hovering over the rug plot across the bottom of each plot includes the metric value at that point as well as the centile that value represents within the data.


# Notes
Input must be a *normalised* VCF. Multiple alts will break the tool and cannot be plotted.

The script currently ignores all metrics except those labelled as FLOAT or INTEGER in the VCF header. This may therefore miss incorrectly labelled annotation metrics (sometimes left as STRING).