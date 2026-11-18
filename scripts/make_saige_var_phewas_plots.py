import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

from adjustText import adjust_text

"""
Variant PheWAS Manhattan plots: one plot per variant_id, phenotypes grouped by category on the x-axis.
Labels phenotypes above the Bonferroni threshold: at most --max_labels_per_category per category and
--max_labels in total, strongest first. Every category without a labeled phenotype also gets its top hit
labeled (not counted toward either cap).

Outputs:
  manhattan_plot_variant_{variant_id}.png - one per variant_id (none if every input is header-only)
  Plotting.log.txt                        - always written
"""

LOG_FILE = 'Plotting.log.txt'


# read in args provided by nextflow channel to process calling python script
def make_arg_parser():
    parser = argparse.ArgumentParser(description="inputs for making single variant plots")

    parser.add_argument("--pheno_descriptions", required=True, help="file containing phenotype name, category and description")
    parser.add_argument("--input_file_list", nargs='+', required=True, help='List of results files')
    parser.add_argument("--max_labels", type=int, default=20, help='Max number of significant phenotypes to label per plot. Default: 20')
    parser.add_argument("--max_labels_per_category", type=int, default=5, help='Max number of significant phenotypes to label per category. Default: 5')
    parser.add_argument("--label_phecode", action='store_true', help='Append the phenotype ID to each label: "DESCRIPTION (phenotype)"')

    return parser


def concatenate_files(files):
    # files are gzipped and tab-separated
    dfs = [pd.read_csv(file, compression='gzip', sep='\t') for file in files]
    return pd.concat(dfs, ignore_index=True)


def log(message):
    print(message)
    with open(LOG_FILE, 'a') as f:
        f.write(message + '\n')


def add_x_positions(df):
    """
    Assign each point an x-position within ±0.4 of its category tick, so points never
    spill into a neighbouring category or past the axis limits.
    Categories are ordered by phenotype count (largest first), like gene-burden v2.
    """
    order = df.groupby('category').size().sort_values(ascending=False, kind='stable').index.tolist()
    df = df.assign(category=pd.Categorical(df['category'], categories=order))
    df = df.sort_values(['category', 'phenotype']).reset_index(drop=True)
    within = df.groupby('category', observed=True).cumcount()
    size = df.groupby('category', observed=True)['phenotype'].transform('size')
    # single-phenotype categories sit on the tick
    offset = np.where(size > 1, within / (size - 1).clip(lower=1) * 0.8 - 0.4, 0.0)
    df['x'] = df['category'].cat.codes + offset
    return df, order


def select_labels(df, thresh, max_labels, max_labels_per_category):
    # significant points, capped per category then in total, plus the top hit of each category without one
    sig = df[df['LOG10P'] > thresh].sort_values('LOG10P', ascending=False)
    sig = sig.groupby('category', observed=True).head(max_labels_per_category).head(max_labels)
    top = df.loc[df.groupby('category', observed=True)['LOG10P'].idxmax()]
    top = top[~top['category'].isin(sig['category'])]
    return pd.concat([sig, top]).drop_duplicates(subset='phenotype')


def plot_variant(df, variant_id, max_labels, max_labels_per_category, label_phecode):
    df, order = add_x_positions(df)
    thresh = -np.log10(0.05 / len(df))
    labeled = select_labels(df, thresh, max_labels, max_labels_per_category)

    palette = dict(zip(order, sns.color_palette('husl', n_colors=len(order))))

    fig, ax = plt.subplots(figsize=(14, 7))
    for category, cat_df in df.groupby('category', observed=True):
        ax.scatter(cat_df['x'], cat_df['LOG10P'], color=palette[category], alpha=0.7, s=10)

    ax.axhline(thresh, color='gray', linestyle='--', linewidth=1)

    texts = []
    for row in labeled.itertuples():
        text = f'{row.description} ({row.phenotype})' if label_phecode else row.description
        texts.append(ax.text(row.x, row.LOG10P, text, color=palette[row.category], size=8,
                             bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=.2', alpha=0.5)))

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=45, ha='right')
    # tick labels in their category's color stand in for a legend
    for tick, category in zip(ax.get_xticklabels(), order):
        tick.set_color(palette[category])
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_ylim(bottom=0, top=max(df['LOG10P'].max(), thresh) * 1.1)
    ax.set_xlabel('Phenotype Category')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'Manhattan Plot for Variant {variant_id}\nBonferroni threshold: p < 0.05 / {len(df)}')

    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

    save_path = f'manhattan_plot_variant_{variant_id}.png'
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return save_path


def main():
    args = make_arg_parser().parse_args()

    # always create the log so the process output exists even when nothing is plotted
    open(LOG_FILE, 'w').close()

    # turn string rep of list into list
    files = [x.strip(',[]') for x in args.input_file_list]
    print(files)
    data = concatenate_files(files)

    data['p_value'] = pd.to_numeric(data['p_value'], errors='coerce')
    data = data.dropna(subset=['p_value'])
    if data.empty:
        log('No results in any input file (e.g. a monomorphic exposure marker dropped by SAIGE); no plots generated.')
        return
    data['LOG10P'] = -np.log10(data['p_value'].clip(lower=1e-300))

    descriptions = pd.read_csv(args.pheno_descriptions)
    descriptions = descriptions.rename(columns={'PHENO': 'phenotype', 'DESCRIPTION': 'description', 'CATEGORY': 'category'})
    # compare IDs as strings so numeric phecodes still match
    descriptions['phenotype'] = descriptions['phenotype'].astype(str)
    data['phenotype'] = data['phenotype'].astype(str)
    data = data.merge(descriptions[['phenotype', 'description', 'category']], on='phenotype', how='left')
    # phenotypes missing from the descriptions file get their own category and are labeled by ID
    data['category'] = data['category'].fillna('Unknown')
    data['description'] = data['description'].fillna(data['phenotype'])

    for variant_id, df in data.groupby('variant_id'):
        save_path = plot_variant(df, variant_id, args.max_labels, args.max_labels_per_category, args.label_phecode)
        log(f'Plot generated for variant_{variant_id} at {save_path}')


if __name__ == '__main__':
    main()
