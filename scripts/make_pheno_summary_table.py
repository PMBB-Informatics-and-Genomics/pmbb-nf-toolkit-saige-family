import pandas as pd
import argparse as ap

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    # Add non-optional list arguments for phenotypes
    parser.add_argument('-b', '--binPhenotypes', nargs='*', required=False, help='List of phenotypes')
    parser.add_argument('-q', '--quantPhenotypes', nargs='*', required=False, help='List of phenotypes')
    parser.add_argument('-t', '--survivalPhenotypes', nargs='*', required=False, help='List of phenotypes')

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

    # Covariates used in SAIGE step 1: samples missing any of these are dropped by SAIGE,
    # so they are excluded here too so that counts match the analyzed sample set
    parser.add_argument('--covars', nargs='*', default=[], help='Covariates (cat + cont) for non sex-stratified cohorts')
    parser.add_argument('--sexStratCovars', nargs='*', default=[], help='Covariates (cat + cont) for sex-stratified cohorts')
    parser.add_argument('--sexStratCohorts', nargs='*', default=[], help='List of sex-stratified cohorts')
    parser.add_argument('--eventTimeCol', required=False, default=None, help='Event time column (survival phenotypes only)')

    return parser

def clean_list(values):
    # Nextflow passes '[]' for empty lists
    return [v for v in (values or []) if v not in ('[]', '')]

args = make_arg_parser().parse_args()

cohorts = args.cohorts
bin_phenos = clean_list(args.binPhenotypes)
quant_phenos = clean_list(args.quantPhenotypes)
survival_phenos = clean_list(args.survivalPhenotypes)
output_dir = args.outDir
id_col = args.id
step1_fam = args.step1Fam
step2_fam = args.step2Fam
step2bgen_sample = args.step2bgen_sample
covars = clean_list(args.covars)
sex_strat_covars = clean_list(args.sexStratCovars)
sex_strat_cohorts = clean_list(args.sexStratCohorts)
event_time_col = args.eventTimeCol

all_phenos = [p for p in bin_phenos]
all_phenos.extend(quant_phenos)

if len(all_phenos) == 0:
    print("No binary or quantitative phenotypes indicated.PLEASE MAKE SURE THAT YOU INTEND TO RUN THIS TYPE  OF ANALYSIS (time-to-event)") #checkpoint 1 to check for binary and quant phenotypes

all_phenos.extend(survival_phenos)

if len(all_phenos) == 0:
    raise Exception("No phenotypes (binary,quantitiative,or survival) indicated")  #checkpoint 2 to check for ANY phenotypes


def is_sex_strat_cohort(cohort):
    # mirrors isSexStratCohort() in processes/saige_helpers.nf
    return any(cohort == s or cohort.startswith(f'{s}__') for s in sex_strat_cohorts)

def cohort_covars(cohort):
    return sex_strat_covars if is_sex_strat_cohort(cohort) else covars


step1_fam = pd.read_table(step1_fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str})

if(step2_fam):
    with open(step2_fam) as _f:
        _is_psam = _f.readline().startswith('#')
    if _is_psam:
        _df = pd.read_table(step2_fam, sep='\\s+', dtype=str)
        _df.columns = [c.lstrip('#') for c in _df.columns]
        _iid_col = 'IID' if 'IID' in _df.columns else _df.columns[0]
        step2_samples = pd.Index(_df[_iid_col])
    else:
        _df = pd.read_table(step2_fam, header=None, names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'],
                            sep='\\s+', dtype={'FID': str, 'IID': str})
        step2_samples = _df['IID']
elif(step2bgen_sample):
    step2bgen_sample = pd.read_table(step2bgen_sample, header=0, comment='#', names=['FID', 'IID', 'MISSING','SEX'], index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str}, skiprows=1)
    step2_samples = step2bgen_sample.index
else:
    raise Exception("No Step 2 input file submitted")

data = pd.read_csv(args.data, index_col=id_col, dtype={id_col: str})
sample_table = pd.read_csv(args.samples, index_col=id_col, dtype={id_col: str})

missing_cols = [c for c in set(all_phenos + covars + sex_strat_covars) if c not in data.columns]
if event_time_col and len(survival_phenos) > 0 and event_time_col not in data.columns:
    missing_cols.append(event_time_col)
if missing_cols:
    raise ValueError(f'Columns not found in data_csv: {missing_cols}')

pheno_types = [(p, 'bin') for p in bin_phenos] + [(p, 'quant') for p in quant_phenos] + [(p, 'survival') for p in survival_phenos]

pheno_info = []
exclusion_info = []

for c in cohorts:
    samples = sample_table.index[sample_table[c] == 1]
    pheno_covars = data.loc[data.index.intersection(samples)]

    keep_samples = sorted(set(samples).intersection(pheno_covars.index).intersection(step1_fam.index).intersection(step2_samples))
    cohort_df = pheno_covars.loc[keep_samples]
    use_covars = cohort_covars(c)

    for p, trait_type in pheno_types:
        # SAIGE drops samples missing the phenotype or any covariate (or event time for survival)
        required_covars = list(use_covars)
        if trait_type == 'survival' and event_time_col:
            required_covars.append(event_time_col)

        missing_pheno = cohort_df[p].isna()
        missing_covar = cohort_df[required_covars].isna().any(axis=1) if required_covars else pd.Series(False, index=cohort_df.index)
        missing_covar_only = missing_covar & ~missing_pheno
        analyzed = cohort_df.loc[~(missing_pheno | missing_covar), p]

        row = {'COHORT': c, 'PHENO': p, 'N': int(analyzed.count())}
        if trait_type in ('bin', 'survival'):
            row['Controls'] = int((analyzed == 0).sum())
            row['Cases'] = int((analyzed == 1).sum())
            row['Prevalence'] = analyzed.mean() if len(analyzed) > 0 else float('nan')
        else:
            row.update(analyzed.describe().drop('count').to_dict())
        pheno_info.append(row)

        covar_detail = ';'.join(f'{cv}:{int(cohort_df.loc[~missing_pheno, cv].isna().sum())}'
                                for cv in required_covars
                                if cohort_df.loc[~missing_pheno, cv].isna().any())
        n_total = len(cohort_df)
        n_excluded = int((missing_pheno | missing_covar).sum())
        if row['N'] == 0:
            note = 'No analyzable samples (all samples missing phenotype or covariates)'
        elif n_excluded > 0:
            note = f'Excluded {n_excluded} of {n_total} samples ({100 * n_excluded / n_total:.1f}%) with missing phenotype/covariates'
        else:
            note = ''
        exclusion_info.append({
            'COHORT': c,
            'PHENO': p,
            'TRAIT_TYPE': trait_type,
            'N_COHORT_GENOTYPED': n_total,
            'N_MISSING_PHENO': int(missing_pheno.sum()),
            'N_MISSING_COVAR_ONLY': int(missing_covar_only.sum()),
            'MISSING_COVAR_DETAIL': covar_detail,
            'N_EXCLUDED': n_excluded,
            'N_ANALYZED': row['N'],
            'CASES': row.get('Cases', ''),
            'CONTROLS': row.get('Controls', ''),
            'NOTE': note,
        })
        if note:
            print(f'WARNING [{c}/{p}]: {note}' + (f' (missing covariates: {covar_detail})' if covar_detail else ''))

pheno_info = pd.DataFrame(pheno_info)
# fixed column order: Groovy parsing in SAIGE_PREPROCESSING reads N at index 2 and Cases at index 4
# Controls/Cases always present (NaN for quant phenos); downstream plot scripts use Cases.count() to detect trait type
for col in ['Controls', 'Cases']:
    if col not in pheno_info.columns:
        pheno_info[col] = float('nan')
col_order = ['COHORT', 'PHENO', 'N', 'Controls', 'Cases']
col_order.extend([c for c in pheno_info.columns if c not in col_order])
pheno_info = pheno_info[col_order]

exclusion_info = pd.DataFrame(exclusion_info)

# specify outdir if given
if output_dir:
    outfile = f'{output_dir}/pheno_summaries.csv'
    exclusion_outfile = f'{output_dir}/sample_exclusion_summary.csv'
else:
    outfile = f'pheno_summaries.csv'
    exclusion_outfile = f'sample_exclusion_summary.csv'
pheno_info.to_csv(outfile, index=False)
exclusion_info.to_csv(exclusion_outfile, index=False)
