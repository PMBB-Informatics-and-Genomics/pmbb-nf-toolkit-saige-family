import pandas as pd
import argparse as ap
import os

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")
    
    # Add non-optional list arguments for phenotypes
    parser.add_argument('-b', '--binPhenotypes', nargs='*', required=False, help='List of phenotypes')
    parser.add_argument('-q', '--quantPhenotypes', nargs='*', required=False, help='List of phenotypes')
    
    # Add a non-optional list argument for cohorts
    parser.add_argument('-c', '--cohorts', nargs='+', required=True, help='List of cohorts')
    
    # Add an argument for .fam file used in step 1
    parser.add_argument('--step1Fam', required=True)
    parser.add_argument('--step2Fam', required=False)
    parser.add_argument('--step2bgen_sample', required=False)

    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default=None, help='Path to output directory. Default: current working directory')

    parser.add_argument('-d', '--data', required=True, help='.csv Phenotype and covariate file')
    parser.add_argument('-s', '--samples', required=True, help='.csv of cohort assignments')

    parser.add_argument('-i', '--id', required=True, help='Column with sample IDs')

    return parser

args = make_arg_parser().parse_args()

cohorts = args.cohorts
bin_phenos = args.binPhenotypes
quant_phenos = args.quantPhenotypes
output_dir = args.outDir
id_col = args.id
step1_fam = args.step1Fam
step2_fam = args.step2Fam
step2bgen_sample = args.step2bgen_sample

all_phenos = [p for p in bin_phenos]
all_phenos.extend(quant_phenos)

if len(all_phenos) == 0:
    raise Exception("No binary or quantitative phenotypes indicated")

print(step1_fam)
import os
os.listdir()
os.listdir('Step1')
step1_fam = pd.read_table(step1_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str})

if(step2_fam):
    step2_fam = pd.read_table(step2_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str})
    step2_samples = step2_fam.index
elif(step2bgen_sample):
    step2bgen_sample = pd.read_table(step2bgen_sample, header=0, comment='#', names=['FID', 'IID', 'MISSING','SEX'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str}, skiprows=1)
    step2_samples = step2bgen_sample.index
else:
    raise Exception("No Step 2 input file submitted")
 
pheno_info = []

data = pd.read_csv(args.data, index_col=id_col, dtype={id_col: str})
sample_table = pd.read_csv(args.samples, index_col=id_col, dtype={id_col: str})

for c in cohorts:
    samples = sample_table.index[sample_table[c] == 1]
    pheno_covars = data.loc[data.index.intersection(samples)]

    keep_samples = list(set(samples).intersection(pheno_covars.index).intersection(step1_fam.index).intersection(step2_samples))

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

pheno_info = pd.concat(pheno_info).reset_index(drop=True)
col_order = ['COHORT', 'PHENO', 'N']
col_order.extend([c for c in pheno_info.columns if c not in col_order])
pheno_info = pheno_info[col_order]
# replacing NAs with 0 for Groovy
pheno_info.loc[pheno_info['N'] == 0, ['Controls', 'Cases']] = 0

# specify outdir if given
if output_dir:
    outfile = f'{output_dir}/pheno_summaries.csv'
else:
    outfile = f'pheno_summaries.csv'
pheno_info.to_csv(outfile, index=False)