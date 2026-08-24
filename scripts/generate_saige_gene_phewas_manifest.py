import json
import argparse
from pathlib import Path

import pandas as pd


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description="Generate a JSON manifest of SAIGE Gene Burden PheWAS pipeline outputs."
    )
    parser.add_argument('--params_json', required=True,
                        help='Path to saige_gene_burden_phewas_params.json')
    parser.add_argument('--pheno_summaries', required=True,
                        help='Path to pheno_summaries.csv')
    parser.add_argument('--singles_summary', required=True,
                        help='Path to singles summary CSV')
    parser.add_argument('--regions_summary', required=True,
                        help='Path to regions summary CSV')
    parser.add_argument('--plots_dir', required=True,
                        help='Directory containing PNG plot files')
    parser.add_argument('--plots_manifests_dir', default=None,
                        help='Directory containing per-chromosome plots_manifest CSVs')
    parser.add_argument('--output', default='saige_gene_phewas_manifest.json',
                        help='Output path for the manifest JSON')
    return parser


def to_list(val):
    if isinstance(val, list):
        return val
    if val is None:
        return []
    return [val]


def main():
    args = make_arg_parser().parse_args()

    with open(args.params_json) as f:
        params = json.load(f)

    cpc = params.get('cohorts_phenotypes_chromosomes', params)
    other = params.get('other_parameters', params)

    cohort_list = to_list(cpc.get('cohort_list', []))
    bin_pheno_list = to_list(cpc.get('bin_pheno_list', []))
    quant_pheno_list = to_list(cpc.get('quant_pheno_list', []))
    survival_pheno_list = to_list(cpc.get('survival_pheno_list', []))
    chromosome_list = to_list(cpc.get('chromosome_list', []))
    p_cutoff = other.get('p_cutoff_summarize', 5e-6)

    plots_dir = Path(args.plots_dir)

    # PheWAS regions plots: {cohort}.{chr}.phewas_regions.{plot_type}.png
    regions_plots = {}
    for cohort in cohort_list:
        regions_plots[cohort] = {}
        for chrom in chromosome_list:
            entry = {}
            for png in sorted(plots_dir.glob(f'{cohort}.{chrom}.phewas_regions.*.png')):
                stem = png.stem
                parts = stem.split('.phewas_regions.')
                if len(parts) == 2:
                    plot_type = parts[1]
                    entry[plot_type] = str(png)
            if entry:
                regions_plots[cohort][str(chrom)] = entry

    # Supplement with actual plot paths from per-chromosome manifest CSVs.
    # The glob above matches stub PNGs; real runs produce per-gene PNGs that
    # are catalogued in the per-chr plots_manifest CSVs instead.
    if args.plots_manifests_dir:
        manifests_dir = Path(args.plots_manifests_dir)
        csv_suffix = '.phewas_regions.plots_manifest.csv'
        for csv_file in sorted(manifests_dir.glob(f'*{csv_suffix}')):
            try:
                df = pd.read_csv(csv_file)
            except Exception:
                continue  # empty (stub) or malformed
            if df.empty or 'gene_phewas_plot' not in df.columns or 'cohort' not in df.columns:
                continue
            cohort = str(df['cohort'].iloc[0])
            # Extract chr from filename: strip '{cohort}.' prefix and csv_suffix
            csv_stem = csv_file.name[:-len(csv_suffix)]
            if not csv_stem.startswith(cohort + '.'):
                continue
            chrom = csv_stem[len(cohort) + 1:]
            if cohort not in regions_plots:
                regions_plots[cohort] = {}
            chr_entry = {}
            for plot_path in df['gene_phewas_plot'].dropna().unique():
                plot_name = Path(plot_path).name
                pt_prefix = cohort + '.'
                pt_suffix = '.phewas.png'
                if plot_name.startswith(pt_prefix) and plot_name.endswith(pt_suffix):
                    plot_type = plot_name[len(pt_prefix):-len(pt_suffix)]
                else:
                    plot_type = Path(plot_path).stem
                chr_entry[plot_type] = str(plots_dir / plot_name)
            if chr_entry:
                regions_plots[cohort][chrom] = chr_entry

    manifest = {
        'cohort_list': cohort_list,
        'bin_pheno_list': bin_pheno_list,
        'quant_pheno_list': quant_pheno_list,
        'survival_pheno_list': survival_pheno_list,
        'chromosome_list': [str(c) for c in chromosome_list],
        'p_cutoff_summarize': p_cutoff,
        'singles_summary_csv': args.singles_summary,
        'regions_summary_csv': args.regions_summary,
        'pheno_summaries_csv': args.pheno_summaries,
        'regions_plots': regions_plots,
        'params': params,
    }

    with open(args.output, 'w') as f:
        json.dump(manifest, f, indent=2)

    n_plots = sum(len(v) for v in regions_plots.values())
    print(f"Manifest written to {args.output}")
    print(f"  Cohorts       : {cohort_list}")
    print(f"  Chromosomes   : {chromosome_list}")
    print(f"  Regions plot combos: {n_plots}")


if __name__ == '__main__':
    main()
