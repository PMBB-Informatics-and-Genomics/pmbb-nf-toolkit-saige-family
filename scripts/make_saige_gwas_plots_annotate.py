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
    parser.add_argument('-s', '--sumstats', required=True, help='Path to summary statistics file')
    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default=None,
                        help='Path to output directory. Default: current working directory')
    # Add a non-optional argument for cohort
    parser.add_argument('-ana', '--analysis', required=True,
                        help='which meta-analysis group of cohorts')
    parser.add_argument('-a', '--annot', required=False, default=None)
    parser.add_argument('-t', '--phenoTable', required=True)
    return parser


args = make_arg_parser().parse_args()
cohort = args.cohort
pheno = args.phenotype
annot_file = args.annot
pheno_table = args.phenoTable
sumstats = args.sumstats
output_dir = args.outDir

pheno_df = pd.read_csv(pheno_table)
pheno_df = pheno_df[pheno_df['PHENO'] == pheno]
pheno_df = pheno_df[pheno_df['COHORT'] == cohort]


trait_type = 'bin' if pheno_df['Cases'].count() != 0 else 'quant'

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

# specify outdir if given
if output_dir:
    output_dir += '/' if output_dir[-1] != '/' else ''
    output_manhattan = f'{output_dir}/{cohort}.{pheno}.manhattan_vertical.png'
    output_qq = f'{output_dir}/{cohort}.{pheno}.qq.png'
else:
    output_manhattan = f'{cohort}.{pheno}.manhattan_vertical.png'
    output_qq = f'{cohort}.{pheno}.qq.png'

plot_title = f'SAIGE GWAS Annotated {cohort}: {pheno.replace("_", " ")}'

if trait_type == 'bin':
    plot_title += f'\nCases = {pheno_df.iloc[0].Cases:,.0f}, Controls = {pheno_df.iloc[0].Controls:,.0f}'
else:
    plot_title += f'\nN = {pheno_df.iloc[0].N:,.0f}'

mp = ManhattanPlot(sumstats, test_rows=None, title=plot_title)
mp.load_data()

mp.df['chromosome_noCHR'] = mp.df['chromosome'].astype(str).str.replace('chr', '').astype(int)
mp.df['variant_id'] = mp.df['variant_id'].str.replace('_', ':')

mp.clean_data(col_map={'chromosome_noCHR': '#CHROM', columns_map['POS']: 'POS', columns_map['MarkerID']: 'ID', columns_map['p.value']: 'P'})

print(mp.df)

if annot_file is not None:
    annot_df = pd.read_csv(annot_file)
    annot_df['ID'] = annot_df['Gene']
    annot_df.iloc[:, 0] = annot_df.iloc[:, 0].str.replace('_', ':')
    mp.add_annotations(annot_df, extra_cols=['RSID'])


mp.get_thinned_data()
print(mp.thinned)
print(len(mp.thinned))

annot_thresh = 5E-8 if np.any(mp.thinned['P'].min() < 5E-8) else np.nanquantile(mp.thinned['P'], 10 / len(mp.thinned))

print(annot_thresh)

mp.update_plotting_parameters(vertical=True, sig=annot_thresh if not np.any(mp.thinned['P'] < 5E-8) else 5E-8, 
                              sug=annot_thresh, annot_thresh=annot_thresh, merge_genes=True)

print(mp.df)

mp.full_plot(save=output_manhattan, rep_boost=False, extra_cols={
             'beta':'BETA', 'RSID': 'RSID'}, number_cols=['BETA'], keep_chr_pos=False, with_table_bg = False, with_table_grid= False)

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