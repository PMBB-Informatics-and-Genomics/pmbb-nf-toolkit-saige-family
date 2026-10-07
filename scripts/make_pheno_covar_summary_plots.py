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
    parser.add_argument('-t', '--survivalPhenotypes', nargs='*', required=False, help='List of phenotypes')
    
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

    # pheno_summaries.csv from make_pheno_summary_table.py (complete-case counts)
    parser.add_argument('--phenoSummary', required=True, help='pheno_summaries.csv')

    # Covariates used in SAIGE step 1: samples missing any of these are excluded from the plots
    parser.add_argument('--covars', nargs='*', default=[], help='Covariates (cat + cont) for non sex-stratified cohorts')
    parser.add_argument('--sexStratCovars', nargs='*', default=[], help='Covariates (cat + cont) for sex-stratified cohorts')
    parser.add_argument('--sexStratCohorts', nargs='*', default=[], help='List of sex-stratified cohorts')
    parser.add_argument('--eventTimeCol', required=False, default=None, help='Event time column (unused for plots)')

    return parser

def clean_list(values):
    # Nextflow passes '[]' for empty lists
    return [v for v in (values or []) if v not in ('[]', '')]

args = make_arg_parser().parse_args()

cohorts = args.cohorts
bin_phenos = clean_list(args.binPhenotypes)
quant_phenos = clean_list(args.quantPhenotypes)
survival_phenos = clean_list(args.survivalPhenotypes)
covars = clean_list(args.covars)
sex_strat_covars = clean_list(args.sexStratCovars)
sex_strat_cohorts = clean_list(args.sexStratCohorts)
step1_fam = args.step1Fam
exome_fam = args.exomeFam
output_dir = args.outDir
id_col = args.id

step1_fam = pd.read_table(step1_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str})

if str(exome_fam).endswith('.psam'):
    exome_fam = pd.read_table(exome_fam, sep='\\s+', dtype=str)
    exome_fam.columns = [c.lstrip('#') for c in exome_fam.columns]
    exome_fam = exome_fam.set_index('IID')
else:
    exome_fam = pd.read_table(exome_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str})

def is_sex_strat_cohort(cohort):
    # mirrors isSexStratCohort() in processes/saige_helpers.nf
    return any(cohort == s or cohort.startswith(f'{s}__') for s in sex_strat_cohorts)

def cohort_covars(cohort):
    return sex_strat_covars if is_sex_strat_cohort(cohort) else covars

subDFs = []

data = pd.read_csv(args.data, index_col=id_col, dtype={id_col: str})
sample_table = pd.read_csv(args.samples, index_col=id_col, dtype={id_col: str})

for c in cohorts:
    samples = sample_table.index[sample_table[c] == 1]
    pheno_covars = data.loc[data.index.intersection(samples)]

    keep_samples = list(set(samples).intersection(pheno_covars.index).intersection(step1_fam.index).intersection(exome_fam.index))

    # drop samples missing any covariate (SAIGE excludes them); phenotype NAs are dropped per-plot below
    subDF = pheno_covars.loc[keep_samples].dropna(subset=cohort_covars(c)).copy()
    subDF['COHORT'] = c
    subDFs.append(subDF)

df_for_violinplots = pd.concat(subDFs).sort_values(by='COHORT')

# Case/control counts come from pheno_summaries.csv so they match the analyzed (complete-case) samples
pheno_info = pd.read_csv(args.phenoSummary, dtype={'PHENO': str})
pheno_info = pheno_info[pheno_info['COHORT'].isin(cohorts)].sort_values(by='COHORT')

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
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

for p, subDF in pheno_info.groupby('PHENO'):
    if p in quant_phenos:
        # Bin phenos only for this one
        continue

    print(p)
    subDF = subDF.copy()
    print(subDF)

    if 'Cases' not in subDF.columns:
        subDF['Cases'] = 0
    subDF['Cases'] = subDF['Cases'].fillna(0).astype(float)
    subDF['N'] = subDF['N'].fillna(0).astype(float)
    # cohorts with no analyzable samples (e.g. phenotype all NA) get prevalence 0 instead of dividing by zero
    subDF['Prevalence'] = (subDF['Cases'] / subDF['N'].where(subDF['N'] > 0)).fillna(0)

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
    
    g = sns.barplot(data=subDF, y='Prevalence', x = 'COHORT',
                    hue='COHORT', palette='turbo', ax=axes[1])
    axes[1].set_ylabel('Prevalence')
    axes[1].set_xticks(axes[1].get_xticks())
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=30, ha='right')
    i = 0
    for _, row in subDF.iterrows():
        prev = row['Prevalence']
        label = '{:.2f}%'.format(prev * 100) if row['N'] > 0 else 'N=0'
        axes[1].text(i, prev, label, ha='center', va='bottom')
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