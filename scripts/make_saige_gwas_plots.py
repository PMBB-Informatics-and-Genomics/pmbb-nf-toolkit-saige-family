from manhattan_plot import ManhattanPlot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import os


def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    # Add a non-optional list argument for phenotype
    parser.add_argument('-p', '--phenotype', required=True, help='phenotype')

    # Add a non-optional list argument for cohort
    parser.add_argument('-c', '--cohort', required=True, help='cohort')
    parser.add_argument('-s', '--sumstats', required=True)
    parser.add_argument('-t', '--phenoTable', required=True)
    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default=None,
                        help='Path to output directory. Default: current working directory')
    # parser.add_argument(
    #     '-m', '--mafList', help="comma-separated list of grouptest MAF values given to SAIGE")
    # parser.add_argument('-a', '--annotationList',
    #                     help="comma-separated list of grouptest annotation values given to SAIGE")
    return parser


args = make_arg_parser().parse_args()
cohort = args.cohort
pheno = args.phenotype

sumstats = args.sumstats
pheno_table = args.phenoTable
output_dir = args.outDir
# for exwas
# grouptest_annotation = args.annotationList
# grouptest_maf = args.mafList


pheno_df = pd.read_csv(pheno_table)
pheno_df = pheno_df[pheno_df['PHENO'] == pheno]
pheno_df = pheno_df[pheno_df['COHORT'] == cohort]

trait_type = 'bin' if pheno_df['Cases'].count() != 0 else 'quant'

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
mp.df['chromosome_noCHR'] = mp.df['chromosome'].astype(str).str.replace('chr', '').astype(int)
mp.clean_data(col_map={'chromosome_noCHR': '#CHROM', columns_map['POS']: 'POS', columns_map['MarkerID']: 'ID', columns_map['p.value']: 'P'})

mp.get_thinned_data()
print(mp.thinned)
print(len(mp.thinned))

# edge case protection
if ~np.any(mp.thinned['P'] < 5E-8):
    p_thresh = np.quantile(mp.thinned['P'], 10 / len(mp.thinned))
else:
		p_thresh = 5E-8

mp.update_plotting_parameters(vertical=True, sig=p_thresh, annot_thresh=p_thresh, merge_genes=True)
mp.full_plot(save=output_manhattan, rep_boost=False, extra_cols={
             'beta': 'beta'}, number_cols=['beta'], keep_chr_pos=False)
plt.clf()

print(f"Saved Manhattan plot to: {output_manhattan}")

mp.qq_plot(save=output_qq)

print(f"Saved qq plot to: {output_qq}")

# plotting manifest
def plots_filename(row):
    manhattan_file = f'Plots/{output_manhattan}'
    qq_file = f'Plots/{output_qq}'
    return (manhattan_file, qq_file)

gwas_sumstats_df = pd.read_csv(sumstats,sep='\t')
gwas_sumstats_df.rename(columns=columns_map_inv,inplace=True)
gwas_sumstats_df.insert(0,'cohort',cohort)
gwas_sumstats_df = gwas_sumstats_df.groupby(['cohort','phenotype']).first().reset_index()
gwas_sumstats_df['results_file'] = sumstats
gwas_sumstats_df[['gwas_manhattan', 'gwas_qq']] = gwas_sumstats_df.apply(lambda row: plots_filename(row), axis=1,result_type='expand')
output_file = f"{cohort}.{pheno}.gwas.plots_manifest.csv"
gwas_sumstats_df.to_csv(output_file,index=False)