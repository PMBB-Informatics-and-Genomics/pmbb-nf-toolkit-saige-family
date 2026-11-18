#!/usr/bin/env python
"""
Harmonize one biobank run cohort's SVP exposure-mode sumstats into the long table consumed by
phewas_meta_tests.py.

Inputs:
  --config       JSON written by Nextflow: biobank, run_cohort, meta_cohort, cont_exposure_transform,
                 exposures {canonical: type}, outcomes (list|null), rename_exposures {local: canonical},
                 rename_outcomes {local: canonical}, files {local exposure: basename}, missing {local exposure: path}
  --encoding     exposure_encoding.tsv (exposure, type, transform, scale, n_exposed), keyed by local exposure name
  --sumstats_dir directory holding the files named in config `files`
Outputs: <out_prefix>.harmonized.tsv.gz, <out_prefix>.harmonization_log.csv

Order of operations per exposure: rename -> filter outcomes -> drop exposure == outcome -> drop invalid stats.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

LONG_COLS = ['biobank', 'cohort', 'exposure', 'exposure_type', 'outcome', 'outcome_type', 'beta',
             'standard_error', 'p_value', 'n', 'n_case', 'n_ctrl', 'n_eff']
LOG_COLS = ['biobank', 'run_cohort', 'meta_cohort', 'item_type', 'name', 'action', 'detail']
NUM_COLS = ['beta', 'standard_error', 'p_value', 'n', 'n_case', 'n_ctrl']


class HarmonizationError(Exception):
    pass


def check_encoding(config, encoding, local, rename_exp):
    """Validate the study's declared type / transform against the run cohort's exposure_encoding."""
    canon = rename_exp.get(local, local)
    bb = config['biobank']
    declared = config['exposures'].get(canon)
    if declared is None:
        raise HarmonizationError(f"{bb}: exposure '{canon}' (local '{local}') is not in the study's exposures")
    rows = encoding[encoding['exposure'] == local]
    if rows.empty:
        raise HarmonizationError(f"{bb}: exposure '{local}' has sumstats but no row in exposure_encoding")
    row = rows.iloc[0]
    if row['type'] != declared:
        raise HarmonizationError(
            f"{bb}: exposure '{canon}' declared '{declared}' in study but exposure_encoding says '{row['type']}'")
    want = config['cont_exposure_transform'] if declared == 'continuous' else 'none'
    if row['transform'] != want:
        raise HarmonizationError(
            f"{bb}: exposure '{canon}' cont_exposure_transform mismatch: study expects '{want}' "
            f"but exposure_encoding says '{row['transform']}'")
    return canon, declared


def harmonize(config, encoding, sumstats):
    bb, rc, mc = config['biobank'], config['run_cohort'], config['meta_cohort']
    log = []

    def add(item_type, name, action, detail=''):
        log.append(dict(biobank=bb, run_cohort=rc, meta_cohort=mc, item_type=item_type, name=name,
                        action=action, detail=detail))

    rename_exp = {str(k): str(v) for k, v in (config.get('rename_exposures') or {}).items()}
    rename_out = {str(k): str(v) for k, v in (config.get('rename_outcomes') or {}).items()}
    keep = config.get('outcomes')
    keep = {str(o) for o in keep} if keep is not None else None
    renamed_outcomes = set()
    parts = []

    for local, path in (config.get('missing') or {}).items():
        add('exposure', local, 'missing_file', f'expected {path}')

    for local in config.get('files') or {}:
        if local not in sumstats:
            raise HarmonizationError(f"{bb}: no sumstats supplied for exposure '{local}'")
        df = sumstats[local]
        canon, etype = check_encoding(config, encoding, local, rename_exp)
        if canon != local:
            add('exposure', local, 'renamed', f'{local} -> {canon}')
        if len(df) == 0:
            add('exposure', canon, 'empty_file', 'no rows')
            continue
        found = set(df['exposure'].astype(str))
        if found != {local}:
            raise HarmonizationError(
                f"{bb}: exposure column of file for '{local}' contains {sorted(found)}, expected '{local}'")
        d = df.copy()
        d['outcome'] = d['phenotype'].astype(str)
        for o in set(d['outcome']) & set(rename_out):
            if o not in renamed_outcomes:
                renamed_outcomes.add(o)
                add('outcome', o, 'renamed', f'{o} -> {rename_out[o]}')
        d['outcome'] = d['outcome'].map(lambda o: rename_out.get(o, o))
        if keep is not None:
            m = ~d['outcome'].isin(keep)
            if m.any():
                add('exposure', canon, 'dropped_not_in_study', f'{int(m.sum())} rows')
            d = d[~m]
        m = d['outcome'] == canon
        for o in d.loc[m, 'outcome'].unique():
            add('outcome', o, 'dropped_self_test', f'outcome equals exposure {canon}')
        d = d[~m].copy()
        for c in NUM_COLS:
            d[c] = pd.to_numeric(d[c], errors='coerce') if c in d.columns else np.nan
        bad = d['beta'].isna() | d['standard_error'].isna() | d['p_value'].isna() | ~(d['standard_error'] > 0)
        if bad.any():
            add('exposure', canon, 'dropped_invalid_stats', f'{int(bad.sum())} rows')
        d = d[~bad].copy()
        d['biobank'], d['cohort'] = bb, mc
        d['exposure'], d['exposure_type'] = canon, etype
        parts.append(d)

    if parts:
        long = pd.concat(parts, ignore_index=True)
    else:
        long = pd.DataFrame(columns=LONG_COLS)
    if len(long):
        binary = long['n_case'].notna()
        long['outcome_type'] = np.where(binary, 'binary', 'quantitative')
        long['n_eff'] = np.where(binary, 2.0 / (1.0 / long['n_case'] + 1.0 / long['n_ctrl']), long['n'])
        dup = long.duplicated(['exposure', 'outcome'], keep=False)
        if dup.any():
            ex = long.loc[dup, ['exposure', 'outcome']].drop_duplicates().head(5).values.tolist()
            raise HarmonizationError(f"{bb}: duplicate (exposure, outcome) rows after renaming: {ex}")
    else:
        long['outcome_type'] = pd.Series(dtype=str)
        long['n_eff'] = pd.Series(dtype=float)
    long = long[LONG_COLS].reset_index(drop=True)
    return long, pd.DataFrame(log, columns=LOG_COLS)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', required=True)
    ap.add_argument('--encoding', required=True)
    ap.add_argument('--sumstats_dir', required=True)
    ap.add_argument('--out_prefix', required=True)
    a = ap.parse_args()
    with open(a.config) as fh:
        config = json.load(fh)
    encoding = pd.read_csv(a.encoding, sep='\t', dtype=str)
    try:
        sumstats = {}
        for local, fname in (config.get('files') or {}).items():
            path = os.path.join(a.sumstats_dir, fname)
            if not os.path.exists(path):
                raise HarmonizationError(f"{config['biobank']}: listed file '{fname}' not found in {a.sumstats_dir}")
            sumstats[local] = pd.read_csv(path, sep='\t', dtype=str)
        long, log = harmonize(config, encoding, sumstats)
    except HarmonizationError as e:
        print(f'ERROR: {e}', file=sys.stderr)
        sys.exit(1)
    long.to_csv(f'{a.out_prefix}.harmonized.tsv.gz', sep='\t', index=False, compression='gzip')
    log.to_csv(f'{a.out_prefix}.harmonization_log.csv', index=False)


if __name__ == '__main__':
    main()
