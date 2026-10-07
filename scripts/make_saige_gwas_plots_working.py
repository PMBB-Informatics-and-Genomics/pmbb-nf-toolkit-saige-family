from manhattan_plot import ManhattanPlot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import os


# def make_arg_parser():
#     parser = ap.ArgumentParser(description=".")

#     # Add a non-optional list argument for phenotype
#     parser.add_argument('-p', '--phenotype', required=True, help='phenotype')

#     # Add a non-optional list argument for cohort
#     parser.add_argument('-c', '--cohort', required=True, help='cohort')

#     parser.add_argument('-g', '--geneFile', required=True)
#     parser.add_argument('-s', '--sumstats', required=True)
#     parser.add_argument('-t', '--phenoTable', required=True)
#     # Add an argument for output_directory
#     parser.add_argument('-o', '--outDir', default=None,
#                         help='Path to output directory. Default: current working directory')
#     # parser.add_argument(
#     #     '-m', '--mafList', help="comma-separated list of grouptest MAF values given to SAIGE")
#     # parser.add_argument('-a', '--annotationList',
#     #                     help="comma-separated list of grouptest annotation values given to SAIGE")
#     return parser
wdir = "/path/to/data/0b146b62b1b8f24fce44a4c8a68279"

cohort = 'AFR_ALL'
pheno = 'T2D'
gene_file = "NCBI.gene.loc"
gene_file = "/path/to/data/NCBI38.gene.loc"

sumstats = "T2D.gwas.saige.gz"
pheno_table = "pheno_summaries.csv"
output_dir = None
# for exwas
# grouptest_annotation = args.annotationList
# grouptest_maf = args.mafList


pheno_df = pd.read_csv(pheno_table)
pheno_df = pheno_df[pheno_df['PHENO'] == pheno]
pheno_df = pheno_df[pheno_df['COHORT'] == cohort]
gene_file = pd.read_table(gene_file, header=None, names=[
                          'ID', 'CHR', 'START', 'STOP', 'STRAND', 'GENE'], index_col='GENE')

trait_type = 'bin' if pheno_df['Cases'].count() != 0 else 'quant'
print(pheno_df)
print(trait_type)
print(gene_file)


# specify outdir if given
if output_dir:
    output_manhattan = f'{output_dir}/{cohort}.{pheno}.manhattan_vertical.png'
    output_qq = f'{output_dir}/{cohort}.{pheno}.qq.png'
else:
    output_manhattan = f'{cohort}.{pheno}.manhattan_vertical.png'
    output_qq = f'{cohort}.{pheno}.qq.png'

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
plot_title = f'SAIGE GWAS {cohort}: {pheno.replace("_", " ")}'

if trait_type == 'bin':
    plot_title += f'\nCases = {pheno_df.iloc[0].Cases:,.0f}, Controls = {pheno_df.iloc[0].Controls:,.0f}'
else:
    plot_title += f'\nN = {pheno_df.iloc[0].N:,.0f}'

mp = ManhattanPlot(sumstats, test_rows=None, title=plot_title)
mp.load_data()
mp.df.rename(columns=columns_map_inv, inplace=True)
mp.df.loc[mp.df['POS'] == 'UR', 'POS'] = gene_file.loc[mp.df.loc[mp.df['POS'] == 'UR', 'MarkerID'].str.split(':', expand=True)[0], 'START'].values
mp.df.loc[mp.df['CHR'] == 'UR', 'CHR'] = gene_file.loc[mp.df.loc[mp.df['CHR'] == 'UR', 'MarkerID'].str.split(':', expand=True)[0], 'CHR'].values
mp.clean_data(col_map={'CHR': '#CHROM', 'MarkerID': 'ID', 'p.value': 'P'})
# mp.clean_data(col_map={columns_map['CHR']: '#CHROM', columns_map['POS']: 'POS', columns_map['MarkerID']: 'ID', columns_map['p.value']: 'P'})
mp.get_thinned_data()
mp.thinned = mp.thinned.dropna(subset='P')
print(mp.thinned)
print(len(mp.thinned))
# num_ind_tests = len(mp.df.ID.unique())
# p_thresh = 0.05 / num_ind_tests

# if ~np.any(mp.thinned['P'] < p_thresh):
#     keep_signals = min(10, len(mp.thinned))
#     p_thresh = 10 ** - \
#         np.nanquantile(mp.thinned['ROUNDED_Y'], 1 -
#                        keep_signals / len(mp.thinned))

# print(p_thresh)
# mp.sig_line = p_thresh

if ~np.any(mp.thinned['P'] < 5E-8):
    p_thresh = np.quantile(mp.thinned['P'], 10 / len(mp.thinned))
else:
		p_thresh = 5E-8

# mp.update_plotting_parameters(vertical=True,sig=p_thresh,sug=p_thresh,annot_thresh=p_thresh,merge_genes=True)
mp.update_plotting_parameters(vertical=True, sig=p_thresh, annot_thresh=p_thresh, merge_genes=True)

mp.full_plot(save=output_manhattan, rep_boost=False, extra_cols={
             'BETA': 'BETA'}, number_cols=['BETA'], keep_chr_pos=False)
plt.clf()

print(f"Saved Manhattan plot to: {output_manhattan}")

mp.qq_plot(save=output_qq)

print(f"Saved qq plot to: {output_qq}")
