#!/usr/bin/env python3
"""Forest plots for the top PheWAS-META hits: one row per cohort plus an IVW diamond."""
import argparse
import hashlib
import math
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def safe_name(s):
    return re.sub(r'[^A-Za-z0-9._-]', '_', str(s))


def plot_stem(exposure, outcome, used):
    """'<exposure>.<outcome>' made file-name safe. If that stem is already in `used` (e.g. outcomes 'A/B' and
    'A_B'), a deterministic suffix (first 8 hex digits of the md5 of the raw names) is appended. Adds to `used`."""
    stem = f'{safe_name(exposure)}.{safe_name(outcome)}'
    if stem in used:
        stem += '.' + hashlib.md5(f'{exposure}\t{outcome}'.encode()).hexdigest()[:8]
    used.add(stem)
    return stem


def _fmt(x):
    """2-3 significant figures; NA for missing."""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return 'NA'
    if math.isnan(x):
        return 'NA'
    return f'{x:.2g}'


def forest_title(row):
    """Both meta p-values (either may be meta_p_col, which ranked the hits) plus heterogeneity."""
    return (f"{row['exposure']} → {row['outcome']}\n"
            f"p_ivw={_fmt(row.get('p_ivw'))}  p_stouffer={_fmt(row.get('p_stouffer'))}  "
            f"p_het={_fmt(row.get('p_het'))}  I²={_fmt(row.get('i2_het'))}")


def make_forest_plot(row, effects_row, cohort_order, out_path):
    """row: filtered-file row (exposure, outcome, beta_ivw, se_ivw, p_ivw, p_stouffer, p_het, i2_het);
    effects_row: effects-table row with beta_<c>, se_<c> per cohort. Cohorts with NaN beta are skipped."""
    labels, betas, ses = [], [], []
    for c in cohort_order:
        b, s = effects_row.get(f'beta_{c}', np.nan), effects_row.get(f'se_{c}', np.nan)
        if pd.notna(b):
            labels.append(c)
            betas.append(float(b))
            ses.append(float(s) if pd.notna(s) else np.nan)
    labels.append('IVW meta')
    n = len(labels)
    ys = np.arange(n, 0, -1)

    fig, ax = plt.subplots(figsize=(7, 1.6 + 0.45 * n))
    for y, b, s in zip(ys[:-1], betas, ses):
        # NaN SE: plot the point without a CI whisker
        ax.errorbar(b, y, xerr=None if np.isnan(s) else 1.96 * s, fmt='s', color='#4c72b0', ecolor='#4c72b0', capsize=3, markersize=5)
    bi, si = float(row['beta_ivw']), float(row['se_ivw'])
    ax.fill([bi - 1.96 * si, bi, bi + 1.96 * si, bi], [ys[-1], ys[-1] + 0.3, ys[-1], ys[-1] - 0.3],
            color='#c44e52', zorder=3)
    ax.axvline(0, color='grey', linestyle='--', linewidth=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels)
    ax.set_ylim(0.4, n + 0.6)
    ax.set_xlabel('beta (95% CI)')
    ax.set_title(forest_title(row), fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--analysis', required=True)
    ap.add_argument('--filtered', required=True)
    ap.add_argument('--effects', required=True)
    ap.add_argument('--cohort_order', nargs='+', required=True)
    ap.add_argument('--max_plots', type=int, default=50)
    a = ap.parse_args()

    filt = pd.read_csv(a.filtered, dtype={'exposure': str, 'outcome': str})
    if filt.empty:
        return
    eff = pd.read_csv(a.effects, sep='\t', dtype={'exposure': str, 'outcome': str})
    eff = eff.set_index(['exposure', 'outcome'])
    used = set()
    for _, row in filt.head(a.max_plots).iterrows():
        key = (row['exposure'], row['outcome'])
        if key not in eff.index:
            continue
        stem = plot_stem(row['exposure'], row['outcome'], used)
        make_forest_plot(row, eff.loc[key], a.cohort_order, f"{a.analysis}.forest.{stem}.png")


if __name__ == '__main__':
    main()
