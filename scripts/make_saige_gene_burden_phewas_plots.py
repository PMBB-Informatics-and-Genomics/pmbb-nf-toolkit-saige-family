import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse as ap
from manhattan_plot import ManhattanPlot

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    # Add a non-optional list argument for cohort
    parser.add_argument('-c', '--cohort', required=True, help='cohort')
    parser.add_argument('-s', '--sumstats', required=True)
    parser.add_argument('--cols')
    parser.add_argument('-d', '--descriptions')

    return parser


# Parse script arguments
args = make_arg_parser().parse_args()
sumstats_file = args.sumstats
cohort = args.cohort
description_file = args.descriptions

# Parse column mapping
colnames_file = args.cols
colnames_rows = open(colnames_file).read().splitlines()
print(colnames_rows)
col_map = dict(zip([r.split('=')[0] for r in colnames_rows],
                   [r.split('=')[1] for r in colnames_rows]))
print(col_map)

# Identify column names to be used later
region_col = 'Region' if 'Region' not in col_map.keys() else col_map['Region']
pvalue_col = 'Pvalue' if 'Pvalue' not in col_map.keys() else col_map['Pvalue']
max_maf_col = 'max_MAF' if 'max_MAF' not in col_map.keys() else col_map['max_MAF']
annot_group_col = 'Group' if 'Group' not in col_map.keys() else col_map['Group']

# Read in summary stats file
sumstats = pd.read_table(sumstats_file)
# Log-transform the p values for plotting
sumstats['LOG10P'] = -np.log10(sumstats[pvalue_col])
print(sumstats)

# Add phenotype descriptions to the results file
pheno_df = pd.read_csv(description_file, index_col='PHENO')
pheno_df = pheno_df.reindex(sumstats['phenotype'])
sumstats[['DESCRIPTION', 'CATEGORY']] = pheno_df[['DESCRIPTION', 'CATEGORY']].values
sumstats['DESCRIPTION'] = sumstats['DESCRIPTION'].fillna('NA')
sumstats['CATEGORY'] = sumstats['CATEGORY'].fillna('Not categorized')
sumstats['CATEGORY'] = sumstats['CATEGORY'].replace('Other Condition', 'Not categorized')
sumstats = sumstats.sort_values(by=['CATEGORY', 'PHENO'])

for region, subDF1 in sumstats.groupby(region_col): # Iterate over genes
    for group, subDF2 in subDF1.groupby(annot_group_col): # Iterate over annotation groups
        for maf, subDF3 in subDF2.groupby(max_maf_col): # Iterate over max MAF thresholds
            print(region, group, maf)

            # Make names for output files
            output_file = f'{cohort}.{region}.{group.replace(";","-")}.maf{str(maf).replace(".","-")}.phewas.png'
            output_qq = f'{cohort}.{region}.{group.replace(";","-")}.maf{str(maf).replace(".","-")}.qq.png'
            print(output_file)

            # Set up plt figure and axes
            fig, ax = plt.subplots()
            fig.set_size_inches(10, 6)

            # Add row number column to serve as x-axis value
            subDF3['row_num'] = np.arange(len(subDF3))

            # Get and Plot Bonferroni corrected P Value Threshold
            p_thresh = -np.log10(0.05 / len(subDF3))
            ax.axhline(p_thresh, linestyle='dashed', color='silver')

            # Scatter plot of the data
            # X = index assigned (row_num)
            # Y = log-transformed p-value (LOG10P)
            # HUE = phenotype category from descriptions file input (CATEGORY)
            sns.scatterplot(data=subDF3, x='row_num', y='LOG10P', hue='CATEGORY', 
                            linewidth=0, palette='turbo', ax=ax, legend=False)
            
            # Get location of tick marks for each phenotype category
            category_ticks = subDF3.groupby('CATEGORY')['row_num'].mean()
            ax.set_xticks(category_ticks)
            ax.set_xticklabels(category_ticks.index, rotation=30, ha='right')

            # Cosmetic plot updates (axis limits, title, labels)
            ax.set_ylim(bottom=0)
            ax.set_xlim(-0.5, len(subDF3)-0.5)
            ax.set_xlabel('PheCode Category')

            plot_title = f'Gene-Burden PheWAS for {region} in {cohort}\nMax MAF = {maf}, Annotation Group = {group}'
            ax.set_title(plot_title)

            # Apply tight layout and write to file
            plt.tight_layout()
            plt.savefig(output_file, bbox_inches='tight')
            plt.clf()
            print(subDF3)

            # ManhattanPlot is used here JUST for QQ plotting because it has it built-in
            mp = ManhattanPlot('', title=plot_title)
            mp.df = subDF3.copy().rename(columns={pvalue_col: 'P'})
            mp.qq_plot(save=output_qq)
            plt.clf()


