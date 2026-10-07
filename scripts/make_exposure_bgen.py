import pandas as pd
import numpy as np
import argparse as ap
from statistics import NormalDist

"""
Encode phenotype exposures as dosage "markers" so SAIGE Step 2 can test them
against the Step 1 null models (exposure mode of the variant PheWAS workflow).

Outputs:
  exposures.vcf             - one marker per exposure on pseudo-chromosome 1, FORMAT GT:DS
  exposures.sample          - Oxford-format sample file matching the VCF sample order
  exposure_encoding.tsv     - per-exposure encoding info, including the BETA/SE back-transform scale
  exposure_cohort_sets.csv  - cohort table with one <cohort>__<exposure> column per pair
                              (cohort member AND exposure non-missing)
"""

COHORT_SEP = '__'


def make_arg_parser():
    parser = ap.ArgumentParser(description='Encode phenotype exposures as a dosage VCF for SAIGE Step 2')

    parser.add_argument('-d', '--data', required=True, help='.csv Phenotype and covariate file')
    parser.add_argument('-i', '--id', required=True, help='Column with sample IDs')
    parser.add_argument('-s', '--samples', required=True, help='.csv of cohort assignments ex: PMBB_AFR_ALL')
    parser.add_argument('-c', '--cohorts', nargs='+', required=True, help='Cohorts to create exposure sample sets for')
    parser.add_argument('--step1Fam', required=True, help='Step 1 plink .fam file')
    parser.add_argument('--bin_exposures', nargs='*', default=[], help='Binary (0/1) exposure columns')
    parser.add_argument('--cont_exposures', nargs='*', default=[], help='Continuous exposure columns')
    parser.add_argument('--cont_transform', default='none', choices=['none', 'inverse_normal'],
                        help='Transform applied to continuous exposures before scaling')
    return parser


def inverse_normal_transform(values):
    # rank-based inverse normal transform (Blom offset), applied to non-missing values only
    ranks = values.rank(method='average')
    n = values.notna().sum()
    nd = NormalDist()
    return ranks.map(lambda r: nd.inv_cdf((r - 0.375) / (n + 0.25)) if pd.notna(r) else np.nan)


def encode_binary(name, values):
    observed = set(values.dropna().unique())
    if not observed.issubset({0, 1}):
        raise ValueError(f'Binary exposure {name} has values other than 0/1: {sorted(observed)[:10]}')
    if len(observed) < 2:
        raise ValueError(f'Binary exposure {name} has no variation')
    # exposed = dosage 1 (het) so that BETA is per exposed individual with no rescaling
    info = {'exposure': name, 'type': 'binary', 'transform': 'none',
            'min': 0, 'max': 1, 'scale': 1.0,
            'n_nonmissing': int(values.notna().sum()), 'n_exposed': int((values == 1).sum())}
    return values.astype(float), info


def encode_continuous(name, values, transform):
    if transform == 'inverse_normal':
        values = inverse_normal_transform(values)
    x_min, x_max = values.min(), values.max()
    if not x_max > x_min:
        raise ValueError(f'Continuous exposure {name} has no variation')
    # scale into [0, 2]; BETA_x = BETA_dosage * 2 / (max - min)
    dosage = 2 * (values - x_min) / (x_max - x_min)
    info = {'exposure': name, 'type': 'continuous', 'transform': transform,
            'min': x_min, 'max': x_max, 'scale': 2 / (x_max - x_min),
            'n_nonmissing': int(values.notna().sum()), 'n_exposed': np.nan}
    return dosage, info


def format_genotype(d):
    if pd.isna(d):
        return './.:.'
    gt = '0/0' if d < 0.5 else ('0/1' if d < 1.5 else '1/1')
    return f'{gt}:{d:.6f}'


def write_vcf(dosages, sample_ids, out_file):
    with open(out_file, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write('##contig=<ID=1>\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Encoded exposure dosage">\n')
        f.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + list(sample_ids)) + '\n')
        for pos, name in enumerate(dosages.columns, start=1):
            genotypes = [format_genotype(d) for d in dosages[name].values]
            f.write('\t'.join(['1', str(pos), name, 'A', 'C', '.', 'PASS', '.', 'GT:DS'] + genotypes) + '\n')


def write_sample_file(sample_ids, out_file):
    # Oxford format, matching the layout of the existing BGEN .sample inputs
    with open(out_file, 'w') as f:
        f.write('ID_1 ID_2 missing sex\n')
        f.write('0 0 0 D\n')
        for s in sample_ids:
            f.write(f'0 {s} 0 NA\n')


args = make_arg_parser().parse_args()
id_col = args.id

exposures = args.bin_exposures + args.cont_exposures
if len(exposures) == 0:
    raise ValueError('No exposures given')
if len(set(exposures)) != len(exposures):
    raise ValueError(f'Duplicate exposure names: {exposures}')

data = pd.read_csv(args.data, index_col=id_col, dtype={id_col: str})
samples = pd.read_csv(args.samples, index_col=id_col, dtype={id_col: str})
step1_fam = pd.read_table(args.step1Fam, header=None, comment='#', names=['FID', 'IID', 'MAT', 'PAT', 'SEX', 'PHENO'],
                          index_col='IID', sep='\\s+', dtype={'FID': str, 'IID': str})

missing_cols = [e for e in exposures if e not in data.columns]
if missing_cols:
    raise ValueError(f'Exposures not found in data_csv: {missing_cols}')

keep_samples = data.index.intersection(step1_fam.index)
data = data.loc[keep_samples]

dosage_cols = {}
encoding_rows = []
for name in args.bin_exposures:
    dosage_cols[name], info = encode_binary(name, pd.to_numeric(data[name], errors='coerce'))
    encoding_rows.append(info)
for name in args.cont_exposures:
    dosage_cols[name], info = encode_continuous(name, pd.to_numeric(data[name], errors='coerce'), args.cont_transform)
    encoding_rows.append(info)

dosages = pd.DataFrame(dosage_cols, index=keep_samples)[exposures]

write_vcf(dosages, keep_samples, 'exposures.vcf')
write_sample_file(keep_samples, 'exposures.sample')
encoding = pd.DataFrame(encoding_rows)
encoding['n_exposed'] = encoding['n_exposed'].astype('Int64')
encoding.to_csv('exposure_encoding.tsv', sep='\t', index=False)

# One derived cohort per (cohort, exposure): cohort members with a non-missing exposure
derived = samples.copy()
for cohort in args.cohorts:
    if cohort not in samples.columns:
        raise ValueError(f'Cohort {cohort} not found in cohort_sets')
    in_cohort = samples[cohort] == 1
    for name in exposures:
        has_exposure = samples.index.isin(dosages.index[dosages[name].notna()])
        derived[f'{cohort}{COHORT_SEP}{name}'] = (in_cohort & has_exposure).astype(int)
derived.index.name = id_col
derived.to_csv('exposure_cohort_sets.csv')

print(encoding.to_string(index=False))
