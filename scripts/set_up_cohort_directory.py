import pandas as pd
import argparse as ap

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")
    
    parser.add_argument('-d', '--data', required=True, help='.csv Phenotype and covariate file')
    parser.add_argument('-c', '--cohort', required=True, help='Cohort to set up')
    parser.add_argument('-s', '--samples', required=True, help='.csv of cohort assignments ex: PMBB_AFR_ALL')
    parser.add_argument('-i', '--id', required=True, help='Column with sample IDs')
    parser.add_argument('--step1Fam', required=True)
    parser.add_argument('--step2Fam', required=False)
    parser.add_argument('--step2bgen_sample', required=False)
    return parser

args = make_arg_parser().parse_args()

id_col = args.id
cohort = args.cohort
step1_fam = args.step1Fam
step2_fam = args.step2Fam
step2bgen_sample = args.step2bgen_sample

data = pd.read_csv(args.data, index_col=id_col, dtype={id_col: str})
samples = pd.read_csv(args.samples, index_col=id_col, dtype={id_col: str})

step1_fam = pd.read_table(step1_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str})

if(step2_fam):
    step2_fam = pd.read_table(step2_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str})
    step2_samples = step2_fam.index
elif(step2bgen_sample):
    step2bgen_sample = pd.read_table(step2bgen_sample, header=0, comment='#', names=['FID', 'IID', 'MISSING','SEX'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str}, skiprows=1)
    step2_samples = step2bgen_sample.index
else:
    raise Exception("No Step 2 input file submitted") 

cohort_samples = samples.index[samples[cohort] == 1]
keep_samples = data.index.intersection(samples.index).intersection(step1_fam.index).intersection(step2_samples).intersection(cohort_samples)

data = data.loc[keep_samples]
step1_fam = step1_fam.loc[keep_samples]

data.index.name = id_col
data.to_csv('saige_pheno_covars.txt', sep='\t')

# The FIDs are usually either the IIDs duplicated or all 0
FIDs = step1_fam['FID']
data.insert(0, 'IID', data.index)
if len(FIDs.unique()) == 1:
    data.insert(0, 'FID', 0)
else:
    data.insert(0, 'FID', data.index)

data[['FID', 'IID']].to_csv(f'sample_list.txt', sep=' ', index=False, header=False)