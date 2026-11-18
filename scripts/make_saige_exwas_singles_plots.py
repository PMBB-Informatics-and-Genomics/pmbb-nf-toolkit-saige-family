from manhattan_plot import ManhattanPlot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import sys
import os


def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    # Add a non-optional list argument for phenotype
    parser.add_argument('-p', '--phenotype', required=True, help='phenotype')

    # Add a non-optional list argument for cohort
    parser.add_argument('-c', '--cohort', required=True, help='cohort')

    parser.add_argument('-g', '--geneFile', required=True)
    parser.add_argument('-s', '--sumstats', required=True)
    parser.add_argument('-t', '--phenoTable', required=True)
    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default=None,
                        help='Path to output directory. Default: current working directory')
    parser.add_argument(
        '-m', '--mafList', help="comma-separated list of grouptest MAF values given to SAIGE")
    parser.add_argument('-a', '--annotationList',
                        help="comma-separated list of grouptest annotation values given to SAIGE")
    return parser


args = make_arg_parser().parse_args()
cohort = args.cohort
pheno = args.phenotype
gene_file = args.geneFile

singles_sumstats = args.sumstats
pheno_table = args.phenoTable
output_dir = args.outDir
grouptest_annotation = args.annotationList
grouptest_maf = args.mafList


pheno_df = pd.read_csv(pheno_table)
pheno_df = pheno_df[pheno_df['PHENO'] == pheno]
pheno_df = pheno_df[pheno_df['COHORT'] == cohort]
gene_file = pd.read_table(gene_file, index_col='gene_id', sep='\s+')

trait_type = 'bin' if pheno_df['Cases'].count() != 0 else 'quant'
print(pheno_df)
print(trait_type)
print(gene_file)


# specify outdir if given
if output_dir:
    output_manhattan = f'{output_dir}/{cohort}.{pheno}.singles.manhattan_vertical.png'
    output_qq = f'{output_dir}/{cohort}.{pheno}.singles.qq.png'
else:
    output_manhattan = f'{cohort}.{pheno}.singles.manhattan_vertical.png'
    output_qq = f'{cohort}.{pheno}.singles.qq.png'

# remap columns from column map file
infile = 'colnames.txt'
# Initialize an empty dictionary
columns_map = {}
with open(infile, 'r') as file:
    for line in file:
        # Split the line on '=' and strip whitespace and quotes
        key, value = line.split('=')
        key = key.strip()
        value = value.strip().strip("'")

        # Add to the dictionary
        columns_map[key] = value
# reverse code the columns map
columns_map_inv = {v: k for k, v in columns_map.items()}

# plot_title = f'SAIGE ExWAS Singles {cohort}: {pheno.replace("_", " ")}'
plot_title = f'SAIGE ExWAS Singles {cohort}: {pheno.replace("_", " ")}'

if trait_type == 'bin':
    plot_title += f'\nCases = {pheno_df.iloc[0].Cases:,.0f}, Controls = {pheno_df.iloc[0].Controls:,.0f}'
else:
    plot_title += f'\nN = {pheno_df.iloc[0].N:,.0f}'

mp = ManhattanPlot(singles_sumstats, test_rows=None, title=plot_title)
mp.load_data()
mp.df.rename(columns=columns_map_inv,inplace=True)
# Replace values less than min float precision with min value
mp.df['p.value'] = np.where(mp.df['p.value'].astype(float) < sys.float_info.min, sys.float_info.min, mp.df['p.value'].astype(float))
mp.clean_data(col_map={'CHR': '#CHROM', 'MarkerID': 'ID', 'p.value': 'P'})
# mp.clean_data(col_map={columns_map['CHR']: '#CHROM', columns_map['POS']: 'POS', columns_map['MarkerID']: 'ID', columns_map['p.value']: 'P'})
mp.get_thinned_data()
mp.thinned = mp.thinned.dropna(subset='P')

num_ind_tests = len(mp.df.ID.unique())
if num_ind_tests == 0:
    print('All Single Assoc Tests were Rare')
    open(output_manhattan, 'w+').write('All Single Assoc Tests Were Ultrarare')
    open(output_qq, 'w+').write('All Single Assoc Tests Were Ultrarare')
    open(output_qq.replace('.png', '.csv'), 'w+').write('All Single Assoc Tests Were Ultrarare')
    sys.exit()

p_thresh = 0.05 / num_ind_tests

if ~np.any(mp.thinned['P'] < p_thresh):
    keep_signals = min(10, len(mp.thinned))
    p_thresh = 10 ** -np.nanquantile(mp.thinned['ROUNDED_Y'], 1 - keep_signals / len(mp.thinned))

print(p_thresh)
mp.sig_line = p_thresh

# Replace underscores with colons to avoide confusing the LaTeX formating
mp.thinned['ID'] = mp.thinned['ID'].str.replace('_', ':')

mp.update_plotting_parameters(vertical=False, sig=p_thresh, annot_thresh=p_thresh, merge_genes=True)
mp.full_plot(save=output_manhattan, rep_boost=False, extra_cols={'BETA': 'BETA'}, number_cols=['BETA'], keep_chr_pos=False, with_table_bg = False, with_table_grid= False)
plt.clf()

print(f"Saved Manhattan plot to: {output_manhattan}")

mp.qq_plot(save=output_qq)

print(f"Saved qq plot to: {output_qq}")

# plotting manifest
def plots_filename(row):
    manhattan_file = f'Plots/{output_manhattan}'
    qq_file = f'Plots/{output_qq}'
    return (manhattan_file, qq_file)
        
singles_sumstats_df = pd.read_csv(singles_sumstats,sep='\t')
singles_sumstats_df.rename(columns=columns_map_inv,inplace=True)
singles_sumstats_df.insert(0,'cohort',cohort)
singles_sumstats_df = singles_sumstats_df.groupby(['cohort','phenotype']).first().reset_index()
singles_sumstats_df['results_file'] = singles_sumstats
singles_sumstats_df[['singles_manhattan', 'singles_qq']] = singles_sumstats_df.apply(lambda row: plots_filename(row), axis=1,result_type='expand')
output_file = f"{cohort}.{pheno}.exwas_singles.plots_manifest.csv"
singles_sumstats_df.to_csv(output_file,index=False)