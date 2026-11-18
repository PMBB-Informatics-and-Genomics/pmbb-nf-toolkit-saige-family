import json
import argparse
from pathlib import Path

import pandas as pd


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description="Generate a JSON manifest of PheWAS-META pipeline outputs."
    )
    parser.add_argument('--params_json', required=True, help='Path to saige_phewas_meta_params.json')
    parser.add_argument('--suggestive', required=True, help='Path to phewas_meta_all_suggestive.csv')
    parser.add_argument('--test_counts', required=True, help='Path to phewas_meta_test_counts.csv')
    parser.add_argument('--sample_sizes', required=True, help='Path to phewas_meta_sample_sizes.csv')
    parser.add_argument('--harmonization_log', required=True, help='Path to phewas_meta_harmonization_log.csv')
    parser.add_argument('--cohorts', required=True, help='Path to phewas_meta_cohorts.csv (meta_cohort,biobank,run_cohort)')
    parser.add_argument('--plots_dir', required=True, help='Directory containing PNG plot files')
    parser.add_argument('--output', default='phewas_meta_manifest.json', help='Output path for the manifest JSON')
    return parser


def read_csv_or_empty(path):
    try:
        return pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def cohort_info(cohorts_csv):
    """meta_cohort -> {biobank, run_cohort}, from phewas_meta_cohorts.csv."""
    df = read_csv_or_empty(cohorts_csv)
    if df.empty:
        return {}
    return {r['meta_cohort']: {'biobank': r['biobank'], 'run_cohort': r['run_cohort']} for _, r in df.iterrows()}


def group_plots(plots_dir, analyses, marker):
    """{analysis: [png paths]} for files named <analysis>.<marker>.<...>.png.

    Analysis names may contain '.', so match on the known analysis names (longest first) rather than splitting.
    """
    grouped = {}
    for png in sorted(Path(plots_dir).glob(f'*.{marker}.*.png')):
        for analysis in sorted(analyses, key=len, reverse=True):
            if png.name.startswith(f'{analysis}.{marker}.'):
                grouped.setdefault(analysis, []).append(str(png))
                break
    return grouped


def build_manifest(params_json, suggestive, test_counts, sample_sizes, harmonization_log, plots_dir, cohorts_csv):
    with open(params_json) as f:
        params = json.load(f)
    analyses_param = params.get('analyses') or {}
    info = cohort_info(cohorts_csv)
    # cohorts per analysis from the params; biobank and run cohort from phewas_meta_cohorts.csv (validated to cover every label)
    analyses = {
        str(name): [{'cohort': str(c), 'biobank': info.get(str(c), {}).get('biobank', '?'),
              'run_cohort': info.get(str(c), {}).get('run_cohort', '?')} for c in cohorts]
        for name, cohorts in analyses_param.items()
    }
    return {
        'params': params,
        'analyses': analyses,
        'harmonization_log': harmonization_log,
        'top_hits': suggestive,
        'test_counts': test_counts,
        'sample_sizes': sample_sizes,
        'phewas_plots': group_plots(plots_dir, analyses, 'phewas_plot'),
        'forest_plots': group_plots(plots_dir, analyses, 'forest'),
    }


def main():
    args = make_arg_parser().parse_args()
    manifest = build_manifest(args.params_json, args.suggestive, args.test_counts, args.sample_sizes,
                              args.harmonization_log, args.plots_dir, args.cohorts)
    with open(args.output, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"Manifest written to {args.output}")
    print(f"  Analyses      : {list(manifest['analyses'])}")
    print(f"  PheWAS plots  : {sum(len(v) for v in manifest['phewas_plots'].values())}")
    print(f"  Forest plots  : {sum(len(v) for v in manifest['forest_plots'].values())}")


if __name__ == '__main__':
    main()
