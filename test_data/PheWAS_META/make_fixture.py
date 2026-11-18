#!/usr/bin/env python
"""
Generate the synthetic two-biobank PheWAS-META fixture next to this script (deterministic, seed 20261110).

Layout (mirrors a real SVP exposure-mode launchDir):
  <BB>/Exposure_Genotypes/exposure_encoding.tsv
  <BB>/<cohort>__<local_exposure>/Sumstats/<cohort>__<local_exposure>.chr1.variant_phewas.saige.gz
  <BB>/Summary/saige_variant_phewas_params.json
  BB2_invnorm/   copy of BB2 whose encoding says inverse_normal for exp_cont_local
  expected_meta.tsv, pheno_descriptions.csv

Edge cases embodied:
  BB1 = POP1, POP2 ; BB2 = AFR, EUR
  BB2 names exp_cont 'exp_cont_local' and outcome y_quant1 'yq1_local'      (rename maps)
  exp_bin appears as an outcome in every file; in exp_bin files it is the self-test row (dropped)
  y_extra appears in all BB2 files only                                      (dropped_not_in_study)
  BB2/EUR__exp__dbl is not written                                           (missing_file)
  BB1/POP2__exp__dbl is header-only                                          (empty_file)
  BB1/POP1__exp__dbl, y_bin1 has NA p_value                                  (dropped_invalid_stats)
  BB2_invnorm encoding has transform inverse_normal for exp_cont_local       (transform mismatch vs 'none')

The expected-value computation below is an independent numpy implementation (no repo code imported).
"""
import gzip
import json
import os
import shutil

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))

SVP_EXPOSURE_COLUMNS = [
    'exposure', 'exposure_type', 'exposure_scale', 'phenotype', 'chromosome', 'base_pair_location', 'variant_id',
    'other_allele', 'effect_allele', 'effect_allele_count', 'effect_allele_frequency', 'missing_rate', 'beta',
    'standard_error', 't_statistic', 'variance', 'p_value', 'p_value_na', 'is_spa_test', 'allele_freq_case',
    'allele_freq_ctrl', 'n_case', 'n_ctrl', 'n_case_hom', 'n_case_het', 'n_ctrl_hom', 'n_ctrl_het', 'n']
ENCODING_COLUMNS = ['exposure', 'type', 'transform', 'min', 'max', 'scale', 'n_nonmissing', 'n_exposed']

BIOBANKS = {'BB1': ['POP1', 'POP2'], 'BB2': ['AFR', 'EUR']}
# canonical exposure -> type ; local names per biobank
EXPOSURES = {'exp_bin': 'binary', 'exp_cont': 'continuous', 'exp__dbl': 'binary'}
LOCAL_EXP = {'BB1': {'exp_bin': 'exp_bin', 'exp_cont': 'exp_cont', 'exp__dbl': 'exp__dbl'},
             'BB2': {'exp_bin': 'exp_bin', 'exp_cont': 'exp_cont_local', 'exp__dbl': 'exp__dbl'}}
BIN_OUT = ['y_bin1', 'y_bin2', '008']
QUANT_OUT = ['y_quant1', 'y_quant2']
LOCAL_OUT = {'BB1': {}, 'BB2': {'y_quant1': 'yq1_local'}}  # canonical -> local
EXPECTED = [('exp_bin', 'y_bin1'), ('exp_cont', 'y_quant1'), ('exp__dbl', 'y_bin2')]


def sumstat_rows(rng, bb, exp_local, exp_type, scale):
    outcomes = [(o, 'binary') for o in BIN_OUT] + [(o, 'quant') for o in QUANT_OUT] + [('exp_bin', 'binary')]
    if bb == 'BB2':
        outcomes.append(('y_extra', 'binary'))
    rows = []
    for o, ot in outcomes:
        o_local = LOCAL_OUT[bb].get(o, o)
        se = round(float(rng.uniform(0.05, 0.3)), 6)
        beta = round(float(rng.normal(0.1, 0.2)), 6)
        p = float(2 * norm.sf(abs(beta / se)))  # consistent with the printed beta / se
        r = dict.fromkeys(SVP_EXPOSURE_COLUMNS, 'NA')
        r.update(exposure=exp_local, exposure_type=exp_type, exposure_scale=scale, phenotype=o_local, chromosome=1,
                 base_pair_location=1, variant_id=exp_local, other_allele='A', effect_allele='C',
                 effect_allele_count=100.0, effect_allele_frequency=0.1, missing_rate=0.0,
                 beta=beta, standard_error=se, t_statistic=round(beta / se, 5),
                 variance=round(1 / se ** 2, 3), p_value=p)
        if ot == 'binary':
            nc, nt = int(rng.integers(40, 400)), int(rng.integers(400, 4000))
            r.update(p_value_na=p, is_spa_test='False', allele_freq_case=0.15, allele_freq_ctrl=0.12,
                     n_case=float(nc), n_ctrl=float(nt), n_case_hom=0.0, n_case_het=10.0, n_ctrl_hom=0.0,
                     n_ctrl_het=100.0)
        else:
            r['n'] = float(rng.integers(500, 5000))
        rows.append(r)
    return pd.DataFrame(rows, columns=SVP_EXPOSURE_COLUMNS)


def write_gz(df, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, 'wt', newline='') as fh:
        df.to_csv(fh, sep='\t', index=False, na_rep='NA', lineterminator='\n')


def write_biobank(rng, bb):
    root = os.path.join(HERE, bb)
    enc = []
    for canon, t in EXPOSURES.items():
        enc.append(dict(exposure=LOCAL_EXP[bb][canon], type=t, transform='none', min=0.0 if t == 'binary' else 0.05,
                        max=1.0 if t == 'binary' else 21.1, scale=1.0 if t == 'binary' else 0.0947,
                        n_nonmissing=1000, n_exposed=280 if t == 'binary' else ''))
    os.makedirs(os.path.join(root, 'Exposure_Genotypes'), exist_ok=True)
    pd.DataFrame(enc, columns=ENCODING_COLUMNS).to_csv(
        os.path.join(root, 'Exposure_Genotypes', 'exposure_encoding.tsv'), sep='\t', index=False)
    os.makedirs(os.path.join(root, 'Summary'), exist_ok=True)
    with open(os.path.join(root, 'Summary', 'saige_variant_phewas_params.json'), 'w') as fh:
        json.dump({'cont_exposure_transform': 'none'}, fh, indent=4)
        fh.write('\n')
    for cohort in BIOBANKS[bb]:
        for canon, t in EXPOSURES.items():
            loc = LOCAL_EXP[bb][canon]
            if bb == 'BB2' and cohort == 'EUR' and canon == 'exp__dbl':
                continue  # missing file
            name = f'{cohort}__{loc}'
            df = sumstat_rows(rng, bb, loc, t, 1.0 if t == 'binary' else 0.0947)
            if bb == 'BB1' and cohort == 'POP2' and canon == 'exp__dbl':
                df = df.iloc[0:0]  # header only
            if bb == 'BB1' and cohort == 'POP1' and canon == 'exp__dbl':
                df.loc[df['phenotype'] == 'y_bin1', ['p_value', 'p_value_na']] = np.nan  # NA p_value
            write_gz(df, os.path.join(root, name, 'Sumstats', f'{name}.chr1.variant_phewas.saige.gz'))


def expected_meta():
    """Independent numpy implementation, reading the files that were written."""
    rows = []
    for exp, out in EXPECTED:
        beta, se, p, neff = [], [], [], []
        for bb, cohorts in BIOBANKS.items():
            for cohort in cohorts:
                loc = LOCAL_EXP[bb][exp]
                path = os.path.join(HERE, bb, f'{cohort}__{loc}', 'Sumstats',
                                    f'{cohort}__{loc}.chr1.variant_phewas.saige.gz')
                if not os.path.exists(path):
                    continue
                d = pd.read_csv(path, sep='\t', dtype={'phenotype': str})
                d = d[d['phenotype'] == LOCAL_OUT[bb].get(out, out)]
                if len(d) == 0 or pd.isna(d['p_value'].iloc[0]):
                    continue
                r = d.iloc[0]
                beta.append(r['beta'])
                se.append(r['standard_error'])
                p.append(r['p_value'])
                neff.append(2 / (1 / r['n_case'] + 1 / r['n_ctrl']) if not pd.isna(r['n_case']) else r['n'])
        beta, se, p, neff = map(np.array, (beta, se, p, neff))
        k = len(beta)
        w = 1 / se ** 2
        b_ivw = (w * beta).sum() / w.sum()
        se_ivw = np.sqrt(1 / w.sum())
        z_i = np.sign(beta) * norm.isf(p / 2)
        z = (np.sqrt(neff) * z_i).sum() / np.sqrt(neff.sum())
        q = (w * (beta - b_ivw) ** 2).sum()
        i2 = max(0.0, (q - (k - 1)) / q)
        rows.append(dict(analysis='ALL', exposure=exp, outcome=out, beta_ivw=b_ivw, se_ivw=se_ivw,
                         p_stouffer=2 * norm.cdf(-abs(z)), q_het=q, i2_het=i2))
    return pd.DataFrame(rows)


def main():
    rng = np.random.default_rng(20261110)
    for d in list(BIOBANKS) + ['BB2_invnorm']:
        shutil.rmtree(os.path.join(HERE, d), ignore_errors=True)
    for bb in BIOBANKS:
        write_biobank(rng, bb)
    shutil.copytree(os.path.join(HERE, 'BB2'), os.path.join(HERE, 'BB2_invnorm'))
    enc_path = os.path.join(HERE, 'BB2_invnorm', 'Exposure_Genotypes', 'exposure_encoding.tsv')
    enc = pd.read_csv(enc_path, sep='\t', dtype=str, keep_default_na=False)
    enc.loc[enc['exposure'] == 'exp_cont_local', 'transform'] = 'inverse_normal'
    enc.to_csv(enc_path, sep='\t', index=False)
    pd.DataFrame({'PHENO': BIN_OUT + QUANT_OUT + ['exp_bin'],
                  'DESCRIPTION': ['Binary outcome 1', 'Binary outcome 2', 'Intestinal infectious diseases',
                                  'Quantitative outcome 1', 'Quantitative outcome 2', 'Binary exposure as outcome'],
                  'CATEGORY': ['cat_a', 'cat_a', 'infectious', 'cat_b', 'cat_b', 'cat_c']}).to_csv(
        os.path.join(HERE, 'pheno_descriptions.csv'), index=False)
    expected_meta().to_csv(os.path.join(HERE, 'expected_meta.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
