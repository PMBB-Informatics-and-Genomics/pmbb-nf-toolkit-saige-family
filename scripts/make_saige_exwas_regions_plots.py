from manhattan_plot import ManhattanPlot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import os
import sys


def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    # Add a non-optional list argument for phenotype
    parser.add_argument('-p', '--phenotype', required=True, help='phenotype')

    # Add a non-optional list argument for cohort
    parser.add_argument('-c', '--cohort', required=True, help='cohort')
    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default=None, help='Path to output directory. Default: current working directory')

    parser.add_argument('-g', '--geneFile')
    parser.add_argument('-s', '--sumstats')
    parser.add_argument('-t', '--phenoTable')
    parser.add_argument('-m', '--mafList', help="comma-separated list of grouptest MAF values given to SAIGE")
    parser.add_argument('-a', '--annotationList', help="comma-separated list of grouptest annotation values given to SAIGE")

    return parser


args = make_arg_parser().parse_args()
cohort = args.cohort
pheno = args.phenotype
gene_file = args.geneFile

regions_sumstats = args.sumstats
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

# loop through all grouptest annotations and MAF values
for ann in grouptest_annotation.split(','):
    for maf in grouptest_maf.split(','):
        maf = np.float64(maf)
        # specify outdir if given
        if output_dir:
            output_manhattan = f'{output_dir}/{cohort}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.manhattan.png'
            output_qq = f'{output_dir}/{cohort}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.qq.png'
        else:
            output_manhattan = f'{cohort}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.manhattan.png'
            output_qq = f'{cohort}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.qq.png'

        # plot_title = f'SAIGE ExWAS Regions {cohort}: {pheno.replace("_", " ")}'
        plot_title = f'SAIGE ExWAS Regions {cohort}: {pheno.replace("_", " ")}\nGroup = {ann}, MAF = {maf}'
        if trait_type == 'bin':
            plot_title += f'\nCases = {pheno_df.iloc[0].Cases:,}, Controls = {pheno_df.iloc[0].Controls:,}'
        else:
            plot_title += f'\nN = {pheno_df.iloc[0].N:,}'

        print(f'\n{plot_title}')

        mp = ManhattanPlot(regions_sumstats, test_rows=None, title=plot_title)
        mp.load_data()
        mp.df.rename(columns=columns_map_inv,inplace=True)

        # Replace values less than min float precision with min value
        mp.df['Pvalue_Burden'] = np.where(mp.df['Pvalue_Burden'].astype(float) < sys.float_info.min, sys.float_info.min, mp.df['Pvalue_Burden'].astype(float))
        mp.df = mp.df[(mp.df['Group'] == ann) & (mp.df['max_MAF'] == maf)]
        print(mp.df)
        # not every combination exists. Pass over iteration if DataFrame is empty
        if mp.df.empty:
            print(
                f"The combination of Group = {ann}, MAF = {maf} does not exist")
            continue
        # Add POS And CHR values from gene_file
        mp.df['POS'] = gene_file.reindex(mp.df['Region'])['seq_region_start'].values
        mp.df['CHR'] = gene_file.reindex(mp.df['Region'])['chromosome'].values
        mp.df['SYMBOL'] = gene_file.reindex(mp.df['Region'])['gene_symbol'].values
        # mp.clean_data(col_map={'CHR': '#CHROM', 'BP': 'POS', columns_map['Region']: 'ID', columns_map['Pvalue_Burden']: 'P'})
        mp.clean_data(col_map={'CHR': '#CHROM', 'BP': 'POS', 'SYMBOL': 'ID', 'Pvalue_Burden': 'P'})
        mp.get_thinned_data()
        mp.thinned = mp.thinned.dropna(subset='P')

        num_ind_tests = len(mp.df.ID.unique())
        p_thresh = 0.05 / num_ind_tests

        print(mp.thinned)

        if len(mp.thinned) <= 1:
            continue

        if ~np.any(mp.thinned['P'] < p_thresh):
            keep_signals = min(10, len(mp.thinned))
            p_thresh = 10 ** -np.nanquantile(mp.thinned['ROUNDED_Y'], 1 - keep_signals / len(mp.thinned))

        print(p_thresh)
        mp.sig_line = p_thresh
        mp.update_plotting_parameters(vertical=False, sig=p_thresh, annot_thresh=p_thresh, 
                                      merge_genes=False, ld_block=0)
        mp.full_plot(save=output_manhattan, rep_boost=False, legend_loc='top')
        plt.clf()
        print(f"Saved Manhattan plot to: {output_manhattan}")

        mp.qq_plot(save=output_qq)

        print(f"Saved qq plot to: {output_qq}")
        
# plotting manifest
def plots_filename(row):
    pheno = row['phenotype']
    ann = row['Group']
    maf = row['max_MAF']
    output_manhattan = f'Plots/{cohort}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.manhattan.png'
    output_qq = f'Plots/{cohort}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.qq.png'
    return (output_manhattan, output_qq)
        
regions_sumstats_df = pd.read_csv(regions_sumstats,sep='\t')
regions_sumstats_df.rename(columns=columns_map_inv,inplace=True)
regions_sumstats_df.insert(0,'cohort',cohort)
regions_sumstats_df = regions_sumstats_df.groupby(['cohort','phenotype', 'Group', 'max_MAF']).first().reset_index()
regions_sumstats_df['results_file'] = regions_sumstats
regions_sumstats_df[['regions_manhattan', 'regions_qq']] = regions_sumstats_df.apply(lambda row: plots_filename(row), axis=1,result_type='expand')
output_file = f"{cohort}.{pheno}.exwas_regions.plots_manifest.csv"
regions_sumstats_df.to_csv(output_file,index=False)