import json
import argparse
from pathlib import Path


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description="Generate a JSON manifest of SAIGE Variant PheWAS pipeline outputs."
    )
    parser.add_argument('--params_json', required=True,
                        help='Path to saige_variant_phewas_params.json')
    parser.add_argument('--pheno_summaries', required=True,
                        help='Path to pheno_summaries.csv')
    parser.add_argument('--singles_summary', required=True,
                        help='Path to singles summary CSV')
    parser.add_argument('--plots_dir', required=True,
                        help='Directory containing PNG plot files')
    parser.add_argument('--output', default='saige_variant_phewas_manifest.json',
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

    bin_pheno_list = to_list(cpc.get('bin_pheno_list', []))
    quant_pheno_list = to_list(cpc.get('quant_pheno_list', []))
    survival_pheno_list = to_list(cpc.get('survival_pheno_list', []))
    p_cutoff = other.get('p_cutoff_summarize', 5e-6)

    plots_dir = Path(args.plots_dir)

    # Variant PheWAS plots: {cohort}.manhattan_plot_variant_{variant_id}.png
    variant_plots = {}
    for png in sorted(plots_dir.glob('*.manhattan_plot_variant_*.png')):
        # Extract cohort and variant id: {cohort}.manhattan_plot_variant_{id}.png
        stem = png.stem
        parts = stem.split('.manhattan_plot_variant_', 1)
        if len(parts) == 2:
            cohort_part, variant_id = parts
            variant_key = f'{cohort_part}.{variant_id}'
            variant_plots[variant_key] = str(png)

    manifest = {
        'bin_pheno_list': bin_pheno_list,
        'quant_pheno_list': quant_pheno_list,
        'survival_pheno_list': survival_pheno_list,
        'p_cutoff_summarize': p_cutoff,
        'singles_summary_csv': args.singles_summary,
        'pheno_summaries_csv': args.pheno_summaries,
        'variant_plots': variant_plots,
        'params': params,
    }

    with open(args.output, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"Manifest written to {args.output}")
    print(f"  Variant plots : {len(variant_plots)}")
    if variant_plots:
        print(f"  Variants      : {list(variant_plots.keys())}")


if __name__ == '__main__':
    main()
