import pandas as pd
import numpy as np
import dominate
from dominate.tags import *
import os
import argparse as ap


def make_arg_parser():
    parser = ap.ArgumentParser(description='Generate HTML report for SAIGE GWAS results.')
    parser.add_argument('-p', '--phenotypes', nargs='+', required=True, help='List of phenotypes')
    parser.add_argument('-c', '--cohorts', nargs='+', required=True, help='List of cohorts')
    parser.add_argument('--gwasColnames', required=True,
                        help='File with GWAS column name mappings (KEY=VALUE per line)')
    parser.add_argument('-o', '--outDir', default=None,
                        help='Output directory. Default: current working directory')
    parser.add_argument('--pval', required=True, type=float,
                        help='P-value cutoff for hit tables')
    return parser


def addCssJsToHead():
    link(href='https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/css/bootstrap.min.css',
         rel='stylesheet',
         integrity='sha384-4bw+/aepP/YC94hEpVNVgiZdgIC5+VKNBQNGCHeKRQN+PtmoHDEXuppvnDJzQIu9',
         crossorigin='anonymous')
    script(src='https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/js/bootstrap.bundle.min.js',
           integrity='sha384-HwwvtgBNo3bZJJLYd8oVXjrBZt8cqVSpeBNS5n7C8IVInixGAoxmnlMuBnhbgrkm',
           crossorigin='anonymous')
    meta(charset='utf-8', name='viewport', content='width=device-width')


def makeMainPage(bin_phenos, quant_phenos):
    '''Constructs the main index page with links to individual phenotype report pages.'''
    mainPage = dominate.document(title='GWAS Results Reports')

    with mainPage.head:
        addCssJsToHead()

    if len(bin_phenos) > 0:
        with mainPage:
            h2('Binary Phenotypes:')
            for pheno in bin_phenos:
                a(pheno, href=pheno + '.html')
                br()

    if len(quant_phenos) > 0:
        with mainPage:
            h2('Quantitative Phenotypes:')
            for pheno in quant_phenos:
                a(pheno, href=pheno + '.html')
                br()

    return mainPage


def makeHitsTable(hits):
    '''Renders an HTML hits table, or a no-hits message if the DataFrame is empty.'''
    if len(hits) == 0:
        p(f'No hits met p-value cutoff ({p_cutoff:.2E})', style='color:red;margin-left:10%')
        return

    printHits = hits.copy()
    p(f'P-value cutoff: {p_cutoff:.2E}', style='color:red;margin-left:10%')
    printHits['P'] = printHits['P'].apply(lambda x: f'{x:.2E}')

    with div(cls='table-responsive'):
        with table(cls='table table-striped'):
            with tr():
                for c in hits.columns:
                    th(c, style='border: 2px solid black')
            for _, row in printHits.iterrows():
                with tr():
                    for _, item in row.items():
                        try:
                            td(str(np.round(item, 5)), style='border: 2px solid black')
                        except Exception:
                            td(item, style='border: 2px solid black')


def makeCohortSection(cohort, pheno, is_binary):
    '''
    Generates the results section for one cohort within a phenotype page.
    Includes sample counts/stats, manhattan + QQ plots, and a hits table.
    '''
    # Skip cohort-pheno combos with no output files in src/
    src_files = [f for f in os.listdir('src') if cohort in f and pheno in f]
    if len(src_files) == 0:
        p('No results found for this cohort.', style='color:gray;margin-left:10%')
        return

    # Sample counts / descriptive stats
    try:
        row = counts.loc[(cohort, pheno)]
        if is_binary:
            p(f'Cases: {int(row["Cases"]):,}', style='margin-left:10%')
            p(f'Controls: {int(row["Controls"]):,}', style='margin-left:10%')
        else:
            p(f'N: {int(row["N"]):,}', style='margin-left:10%')
            p(f'Mean: {row["mean"]:.2f}', style='margin-left:10%')
            p(f'SD: {row["std"]:.2f}', style='margin-left:10%')
    except KeyError:
        pass

    # Manhattan and QQ plots
    manhattan_file = f'src/{cohort}.{pheno}.manhattan_vertical.png'
    qq_file = f'src/{cohort}.{pheno}.qq.png'

    if os.path.exists(manhattan_file):
        img(src=manhattan_file, cls='img-fluid')
        br()
    if os.path.exists(qq_file):
        img(src=qq_file, cls='img-fluid')
        br()

    # Hits table from the suggestive summary
    hits_file = 'src/saige_gwas_suggestive.csv'
    if not os.path.exists(hits_file):
        return

    try:
        hits = pd.read_csv(hits_file)
        # Remap actual column names -> canonical SAIGE names, then P-value to 'P'
        hits = hits.rename(columns=gwas_col_map_inv).rename(columns={'p.value': 'P'})
        hits = hits[(hits['cohort'] == cohort) & (hits['phenotype'] == pheno)]
        hits = hits[hits['P'] <= p_cutoff].sort_values('P')
        keep_cols = [c for c in ['MarkerID', 'CHR', 'POS', 'Allele1', 'Allele2',
                                 'AF_Allele2', 'BETA', 'SE', 'P'] if c in hits.columns]
        hits = hits[keep_cols].drop_duplicates()
        makeHitsTable(hits)
    except Exception as e:
        p(f'Could not load hits table: {e}', style='color:orange;margin-left:10%')


# ── Script entry point ────────────────────────────────────────────────────────

args = make_arg_parser().parse_args()
p_cutoff = args.pval
cohorts = args.cohorts
all_phenos = args.phenotypes
output_dir = args.outDir

# Parse GWAS column name mapping: KEY=actual_col_name -> {actual: canonical}
gwas_col_map = {}
with open(args.gwasColnames) as fh:
    for line in fh:
        line = line.strip()
        if '=' not in line:
            continue
        key, val = line.split('=', 1)
        gwas_col_map[key.strip()] = val.strip().strip("'")
gwas_col_map_inv = {v: k for k, v in gwas_col_map.items()}

# Load phenotype summary table (COHORT, PHENO indexed)
counts = pd.read_csv('src/pheno_summaries.csv', dtype={'PHENO': str})
counts = counts.set_index(['COHORT', 'PHENO'])

# Derive binary and quant pheno lists from the summary table
bin_phenos = (counts.dropna(subset=['Cases']).index.get_level_values('PHENO').unique().tolist()
              if 'Cases' in counts.columns else [])
quant_phenos = (counts.dropna(subset=['mean']).index.get_level_values('PHENO').unique().tolist()
                if 'mean' in counts.columns else [])

# Use cohorts from the table (respects filtering done in preprocessing)
cohorts = sorted(counts.index.get_level_values('COHORT').unique().tolist()) or cohorts

print(f'Cohorts: {cohorts}')
print(f'Binary phenotypes: {bin_phenos}')
print(f'Quantitative phenotypes: {quant_phenos}')
print(f'P-value cutoff: {p_cutoff:.2E}')

# Generate index page
index_doc = makeMainPage(bin_phenos, quant_phenos)
out_index = f'{output_dir}/index.html' if output_dir else 'index.html'
open(out_index, 'w').write(str(index_doc))

# Generate one page per phenotype
for pheno in list(bin_phenos) + list(quant_phenos):
    is_binary = pheno in bin_phenos

    codePage = dominate.document(title=f'{pheno}: GWAS Results Report')
    with codePage.head:
        addCssJsToHead()

    with codePage:
        with div(cls='container'):
            header(a('Back to Index', href='index.html'))
            h1(pheno)

            # Phenotype-level summary plot
            summary_plot = (f'src/{pheno}.barplots.png' if is_binary
                            else f'src/{pheno}.violinplot.png')
            if os.path.exists(summary_plot):
                img(src=summary_plot, cls='img-fluid')

            # Per-cohort collapsible sections
            with div():
                with p():
                    for cohort in cohorts:
                        butt = button(f'Cohort: {cohort}', cls='btn btn-light')
                        butt['type'] = 'button'
                        butt['data-bs-toggle'] = 'collapse'
                        butt['data-bs-target'] = f'#collapse{cohort}'
                        butt['aria-expanded'] = 'false'
                        butt['aria-controls'] = f'collapse{cohort}'

                for cohort in cohorts:
                    with div(id=f'collapse{cohort}', cls='collapse multi-collapse'):
                        with div(cls='card card-body'):
                            h3(f'Cohort: {cohort}')
                            makeCohortSection(cohort, pheno, is_binary)

    out_pheno = f'{output_dir}/{pheno}.html' if output_dir else f'{pheno}.html'
    open(out_pheno, 'w').write(str(codePage))
