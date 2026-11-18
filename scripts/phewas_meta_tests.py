#!/usr/bin/env python
"""
Meta-analysis of exposure-mode PheWAS results across cohorts.

Input: long table (one row per exposure x outcome x cohort) with columns
  biobank, cohort, exposure, exposure_type, outcome, outcome_type, beta, standard_error, p_value,
  n, n_case, n_ctrl, n_eff
Per key (default exposure x outcome) computes fixed-effect inverse-variance meta-analysis, sample-size
weighted signed Stouffer Z, and Cochran's Q / I^2 heterogeneity.

Outputs (prefix = {analysis}.phewas_meta):
  .tsv.gz              - one row per key, sorted by exposure then --p_col
  .filtered.csv        - rows with p_col <= pthresh and N_studies >= min_studies, `analysis` column first
  .effects.tsv.gz      - per-cohort beta/se/p/n_eff, wide
  .test_counts.csv     - analysis,exposure,n_tests
  .sample_sizes.csv    - analysis,exposure,outcome_type,max_n_eff,max_N_studies
  .plot_input.tsv.gz   - variant_id,phenotype,p_value for make_saige_var_phewas_plots.py
"""

import argparse

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm


META_COLS = ['N_studies', 'direction', 'n_total', 'n_case_total', 'n_ctrl_total', 'n_eff_total',
             'beta_ivw', 'se_ivw', 'p_ivw', 'z_stouffer', 'p_stouffer', 'q_het', 'p_het', 'i2_het']
ID_COLS = ['exposure', 'outcome', 'cohort', 'biobank']
STUDY_COLS = ['beta', 'standard_error', 'p_value', 'n', 'n_case', 'n_ctrl', 'n_eff']


def meta_analyze(long_df, key_cols, cohort_order):
    """Returns (meta_df, effects_df); see module docstring.

    Vectorized: rows are numbered by key (sorted, as groupby(sort=True)) and every per-key sum is a bincount.
    Rows with missing beta/SE or SE <= 0 are dropped first; a key x cohort pair may appear only once."""
    meta_cols = key_cols + ['exposure_type', 'outcome_type'] + META_COLS
    eff_cols = key_cols + [f'{s}_{c}' for c in cohort_order for s in ('beta', 'se', 'p', 'n_eff')]

    df = long_df.dropna(subset=['beta', 'standard_error'])
    df = df[df['standard_error'] > 0]
    # rows with a missing key value are skipped, as groupby(dropna=True) does
    df = df.dropna(subset=key_cols).reset_index(drop=True)
    if df.empty:
        return pd.DataFrame(columns=meta_cols), pd.DataFrame(columns=eff_cols)

    codes = df.groupby(key_cols, sort=True).ngroup().to_numpy()
    n_keys = int(codes.max()) + 1
    cohort_codes, cohort_uniques = pd.factorize(df['cohort'], use_na_sentinel=False)
    dup_mask = pd.Series(codes * len(cohort_uniques) + cohort_codes).duplicated(keep=False).to_numpy()
    if dup_mask.any():
        dup_cols = key_cols + ['cohort']
        dups = df.loc[dup_mask, dup_cols].drop_duplicates()
        shown = [tuple(r) for r in dups.head(5).itertuples(index=False, name=None)]
        raise ValueError(f'{len(dups)} duplicate ({", ".join(dup_cols)}) combinations in input '
                         f'(harmonization bug or shared cohort label?); first {len(shown)}: {shown}')

    def ksum(x):
        return np.bincount(codes, weights=x, minlength=n_keys)

    def ksum_or_nan(col):
        # sum of the non-missing values; NaN when a key has none
        x = df[col].to_numpy(float)
        has = ~np.isnan(x)
        total = ksum(np.where(has, x, 0.0))
        total[ksum(has.astype(float)) == 0] = np.nan
        return total

    beta = df['beta'].to_numpy(float)
    se = df['standard_error'].to_numpy(float)
    p = df['p_value'].to_numpy(float)
    n_eff = df['n_eff'].to_numpy(float)
    k = np.bincount(codes, minlength=n_keys)

    # fixed-effect inverse-variance
    w = 1.0 / se ** 2
    sw = ksum(w)
    beta_ivw = ksum(w * beta) / sw
    se_ivw = np.sqrt(1.0 / sw)
    p_ivw = 2 * norm.sf(np.abs(beta_ivw / se_ivw))

    # signed, sample-size weighted Stouffer over rows with finite Z and finite n_eff > 0;
    # isf keeps p == 0 finite (clipped at 1e-300)
    z_i = np.sign(beta) * norm.isf(np.clip(p, 1e-300, 1.0) / 2)
    with np.errstate(invalid='ignore'):
        ok = np.isfinite(z_i) & np.isfinite(n_eff) & (n_eff > 0)
    n_eff_ok = np.where(ok, n_eff, 0.0)
    s_wz = ksum(np.sqrt(n_eff_ok) * np.where(ok, z_i, 0.0))
    s_n = ksum(n_eff_ok)
    has_z = ksum(ok.astype(float)) > 0
    z = np.full(n_keys, np.nan)
    z[has_z] = s_wz[has_z] / np.sqrt(s_n[has_z])
    p_st = 2 * norm.sf(np.abs(z))

    # Cochran's Q (deviations from each key's IVW estimate) and I^2 (floored at 0); NaN for a single study
    q = ksum(w * (beta - beta_ivw[codes]) ** 2)
    multi = k >= 2
    p_het = np.full(n_keys, np.nan)
    p_het[multi] = chi2.sf(q[multi], k[multi] - 1)
    i2 = np.full(n_keys, np.nan)
    pos = multi & (q > 0)
    i2[pos] = np.maximum(0.0, (q[pos] - (k[pos] - 1)) / q[pos])
    i2[multi & ~(q > 0)] = 0.0
    q = np.where(multi, q, np.nan)

    # per-cohort cells in cohort_order (a cohort not in cohort_order counts in the stats but gets no column)
    n_c = len(cohort_order)
    col = pd.Index(list(cohort_order)).get_indexer(df['cohort'])
    in_order = col >= 0
    rows, cols = codes[in_order], col[in_order]
    sign = np.full((n_keys, max(n_c, 1)), '?', dtype='<U1')
    b_in = beta[in_order]
    sign[rows, cols] = np.where(b_in > 0, '+', np.where(b_in < 0, '-', '0'))
    direction = sign.view(f'<U{n_c}').ravel().astype(object) if n_c else np.full(n_keys, '', dtype=object)

    # one row per key, in key order; exposure_type / outcome_type come from the key's first input row
    _, first = np.unique(codes, return_index=True)
    keys = df.iloc[first].reset_index(drop=True)
    meta = keys[key_cols + ['exposure_type', 'outcome_type']].copy()
    stats = dict(N_studies=k, direction=direction, n_total=ksum_or_nan('n'), n_case_total=ksum_or_nan('n_case'),
                 n_ctrl_total=ksum_or_nan('n_ctrl'), n_eff_total=ksum_or_nan('n_eff'), beta_ivw=beta_ivw,
                 se_ivw=se_ivw, p_ivw=p_ivw, z_stouffer=z, p_stouffer=p_st, q_het=q, p_het=p_het, i2_het=i2)
    meta = pd.concat([meta, pd.DataFrame(stats)], axis=1)

    wide = {}
    for j, c in enumerate(cohort_order):
        for stat, x in (('beta', beta), ('se', se), ('p', p), ('n_eff', n_eff)):
            a = np.full(n_keys, np.nan)
            sel = cols == j
            a[rows[sel]] = x[in_order][sel]
            wide[f'{stat}_{c}'] = a
    eff = pd.concat([keys[key_cols], pd.DataFrame(wide)], axis=1)
    return meta[meta_cols], eff[eff_cols]


def make_arg_parser():
    parser = argparse.ArgumentParser(description='Meta-analyze exposure-mode PheWAS results across cohorts')
    parser.add_argument('--analysis', required=True, help='analysis name, used as output prefix')
    parser.add_argument('--cohort_order', nargs='+', required=True, help='meta cohort labels; order of direction string and effects columns')
    parser.add_argument('--inputs', nargs='+', required=True, help='harmonized long tables (tsv.gz)')
    parser.add_argument('--pthresh', type=float, required=True, help='p threshold for the filtered output')
    parser.add_argument('--min_studies', type=int, default=2, help='min N_studies for the filtered output. Default: 2')
    parser.add_argument('--p_col', default='p_ivw', help='p-value column used for sorting/filtering/plots. Default: p_ivw')
    parser.add_argument('--key_cols', nargs='+', default=['exposure', 'outcome'], help='columns defining a test. Default: exposure outcome')
    return parser


def main():
    args = make_arg_parser().parse_args()
    dtypes = {c: str for c in ID_COLS}
    long_df = pd.concat([pd.read_csv(f, sep='\t', dtype=dtypes) for f in args.inputs], ignore_index=True)

    meta_df, effects_df = meta_analyze(long_df, args.key_cols, args.cohort_order)
    meta_df = meta_df.sort_values(['exposure', args.p_col], kind='mergesort').reset_index(drop=True)
    pre = f'{args.analysis}.phewas_meta'

    meta_df.to_csv(f'{pre}.tsv.gz', sep='\t', index=False)
    effects_df.to_csv(f'{pre}.effects.tsv.gz', sep='\t', index=False)

    filt = meta_df[(meta_df[args.p_col] <= args.pthresh) & (meta_df['N_studies'] >= args.min_studies)]
    filt = filt.sort_values(args.p_col, kind='mergesort')
    filt.insert(0, 'analysis', args.analysis)
    filt.to_csv(f'{pre}.filtered.csv', index=False)

    counts = meta_df.groupby('exposure').size().rename('n_tests').reset_index()
    counts.insert(0, 'analysis', args.analysis)
    counts.to_csv(f'{pre}.test_counts.csv', index=False)

    sizes = (meta_df.groupby(['exposure', 'outcome_type'])
             .agg(max_n_eff=('n_eff_total', 'max'), max_N_studies=('N_studies', 'max')).reset_index())
    sizes.insert(0, 'analysis', args.analysis)
    sizes.to_csv(f'{pre}.sample_sizes.csv', index=False)

    plot = pd.DataFrame({'variant_id': meta_df['exposure'], 'phenotype': meta_df['outcome'],
                         'p_value': meta_df[args.p_col]})
    plot.to_csv(f'{pre}.plot_input.tsv.gz', sep='\t', index=False)


if __name__ == '__main__':
    main()
