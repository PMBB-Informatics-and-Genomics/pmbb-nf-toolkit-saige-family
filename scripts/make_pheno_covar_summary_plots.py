import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse as ap
import os

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")
    
    # Add non-optional list arguments for phenotypes
    parser.add_argument('-b', '--binPhenotypes', nargs='*', required=False, help='List of phenotypes')
    parser.add_argument('-q', '--quantPhenotypes', nargs='*', required=False, help='List of phenotypes')
    
    # Add a non-optional list argument for cohorts
    parser.add_argument('-c', '--cohorts', nargs='+', required=True, help='List of cohorts')

    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default=None, help='Path to output directory. Default: current working directory')
    
    # Add an argument for .fam file used in step 1
    parser.add_argument('--step1Fam', required=True)

    # Add an argument for .fam exome file
    parser.add_argument('--exomeFam', required=True)

    parser.add_argument('-d', '--data', required=True, help='.csv Phenotype and covariate file')
    parser.add_argument('-s', '--samples', required=True, help='.csv of cohort assignments')

    parser.add_argument('-i', '--id', required=True, help='Column with sample IDs')

    return parser

args = make_arg_parser().parse_args()

cohorts = args.cohorts
bin_phenos = args.binPhenotypes
quant_phenos = args.quantPhenotypes
step1_fam = args.step1Fam
exome_fam = args.exomeFam
output_dir = args.outDir
id_col = args.id

step1_fam = pd.read_table(step1_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str})
exome_fam = pd.read_table(exome_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str})

pheno_info = []

subDFs = []

data = pd.read_csv(args.data, index_col=id_col, dtype={id_col: str})
sample_table = pd.read_csv(args.samples, index_col=id_col, dtype={id_col: str})

for c in cohorts:
    samples = sample_table.index[sample_table[c] == 1]
    pheno_covars = data.loc[data.index.intersection(samples)]

    keep_samples = list(set(samples).intersection(pheno_covars.index).intersection(step1_fam.index).intersection(exome_fam.index))

    if len(bin_phenos) > 0:
        bin_pheno_info = pheno_covars.loc[keep_samples, bin_phenos].apply(lambda x: x.value_counts(), result_type='expand').transpose().rename(columns={0: 'Controls', 1: 'Cases'})
        bin_pheno_info['N'] = pheno_covars.loc[keep_samples, bin_phenos].count()
        bin_pheno_info['Prevalence'] = pheno_covars.loc[keep_samples, bin_phenos].mean()
        bin_pheno_info.index.name = 'PHENO'
        bin_pheno_info = bin_pheno_info.reset_index()
        bin_pheno_info['COHORT'] = c
        pheno_info.append(bin_pheno_info)

    if len(quant_phenos) > 0:
        quant_pheno_info = pheno_covars.loc[keep_samples, quant_phenos].describe().drop('count').transpose()
        quant_pheno_info['N'] = pheno_covars.loc[keep_samples, quant_phenos].count()
        quant_pheno_info.index.name = 'PHENO'
        quant_pheno_info = quant_pheno_info.reset_index()
        quant_pheno_info['COHORT'] = c
        pheno_info.append(quant_pheno_info)

    subDF = pheno_covars.loc[keep_samples].copy()
    subDF['COHORT'] = c
    subDFs.append(subDF)

df_for_violinplots = pd.concat(subDFs).sort_values(by='COHORT')
pheno_info = pd.concat(pheno_info).reset_index(drop=True).sort_values(by='COHORT')

for p in quant_phenos:
    print(df_for_violinplots[[p, 'COHORT']])
    sns.violinplot(data=df_for_violinplots.dropna(subset=[p]), y=p, x='COHORT', palette='turbo', hue='COHORT')
    plt.gca().set_xticks(plt.gca().get_xticks())
    plt.gca().set_xticklabels(plt.gca().get_xticklabels(), rotation=30, ha='right')
    # specify output directory if given
    if output_dir:
        outfile=f'{output_dir}/{p}.violinplot.png'
    else:
        outfile=f'{p}.violinplot.png'
    plt.savefig(outfile)
    plt.clf()

for p, subDF in pheno_info.groupby('PHENO'):
    if p in quant_phenos:
        # Bin phenos only for this one
        continue

    print(p)
    print(subDF)

    subDF['Cases'] = subDF['Cases'].fillna(0)

    fig, axes = plt.subplots(ncols=2)
    fig.set_size_inches(10, 5)

    g = sns.barplot(data=subDF, y='Cases', x='COHORT', 
                    hue='COHORT', palette='turbo', ax=axes[0])
    g.set_yscale("log")
    axes[0].set_xticks(axes[0].get_xticks())
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=30, ha='right')
    i = 0
    for _, row in subDF.iterrows():
        axes[0].text(i, row['Cases'], '{:,}'.format(int(row['Cases'])), ha='center', va='bottom')
        i += 1
    
    g = sns.barplot(data=subDF, y=subDF['Cases'] / subDF['N'], x = 'COHORT',
                    hue='COHORT', palette='turbo', ax=axes[1])
    axes[1].set_ylabel('Prevalence')
    axes[1].set_xticks(axes[1].get_xticks())
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=30, ha='right')
    i = 0
    for _, row in subDF.iterrows():
        prev = row['Cases'] / row['N']
        axes[1].text(i, prev, '{:.2f}%'.format(prev * 100), ha='center', va='bottom')
        i += 1
    plt.tight_layout()
    # specify outdir if given
    if output_dir:
        outfile=f'{output_dir}/{p}.barplots.png'
    else:
        outfile=f'{p}.barplots.png'
    plt.savefig(outfile)
    # plt.savefig(f'{output_dir}/{p}.barplots.png')
    plt.clf()