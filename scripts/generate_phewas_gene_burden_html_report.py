import pandas as pd
import numpy as np
import dominate
from dominate.tags import *
import os
import sys
import argparse as ap

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")
    
    # Add a non-optional list argument for phenotypes
    parser.add_argument('-p', '--phenotypes', nargs='+', required=True, help='List of phenotypes')
    
    # Add a non-optional list argument for cohorts
    parser.add_argument('-c', '--cohorts', nargs='+', required=True, help='List of cohorts')

    parser.add_argument('--regionColnames')
    parser.add_argument('--singlesColnames')

    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default=None, help='Path to output directory. Default: current working directory')
    
    # Add an argument for p-value cutoff
    parser.add_argument('--pval', required=True, help='P-Value Cutoff for Hit Tables', type=float)
    
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
    mainPage = dominate.document(title = 'PheWAS Results Reports')

    with mainPage.head:
        addCssJsToHead()
    
    with mainPage:
        a('Methods:', href='src/saige_exwas_methods.html')

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
    if (len(hits) == 0):
        p('No Hits Met User-Specified P-Value Filter ({p:.2E})'.format(p=p_cutoff), style='color:red;margin-left:10%')
        return
    printHits = hits.copy()
    p('User-Specified P-Value Filter: {p:.2E}'.format(p=p_cutoff), style='color:red;margin-left:10%')
    printHits['P'] = printHits['P'].apply(lambda x: '{:.2E}'.format(x))
    with div(cls='table-responsive'):
        with table(cls="table table-striped"):
            with tr():
                for c in hits.columns:
                    th(c, style='border: 2px solid black')
            for _, row in printHits.iterrows():
                with tr():
                    for _, item in row.items():
                        try:
                            td(str(np.round(item, 5)), style='border: 2px solid black')
                        except:
                            td(item, style='border: 2px solid black')


def makeSingleBinaryReport(cohort, pheno):
    # Get Source Files (hits table and plots)
    src_files = ['src/' + f for f in os.listdir('src') if cohort in f and pheno in f]
    if len(src_files) < 2:
        return
    if counts.loc[(cohort, pheno), 'Cases'] < 50:
        return
    
    p('Cases: ' + str(counts.loc[(cohort, pheno), 'Cases'].astype(int)), style='margin-left:10%')
    p('Controls: ' + str(counts.loc[(cohort, pheno), 'Controls'].astype(int)), style='margin-left:10%')

    region_html_id = 'collapse' + '-'.join([cohort, 'regions'])
    singles_html_id = 'collapse' + '-'.join([cohort, 'singles'])

    with p():
        butt = button('SAIGE-GENE Region Results (' + cohort + ')', cls='btn btn-light')
        butt['type'] = 'button'
        butt['data-bs-toggle'] = 'collapse'
        butt['data-bs-target'] = '#' + region_html_id
        butt['aria-expanded'] = 'false'
        butt['aria-controls'] = region_html_id
        
        butt = button('SAIGE-GENE Singles Results (' + cohort + ')', cls='btn btn-light')
        butt['type'] = 'button'
        butt['data-bs-toggle'] = 'collapse'
        butt['data-bs-target'] = '#' + singles_html_id
        butt['aria-expanded'] = 'false'
        butt['aria-controls'] = singles_html_id

    with div(id=region_html_id, cls='collapse multi-collapse'):
        with div(cls='card card-body'):
            h3('SAIGE-GENE Regions:')
            manhattan_files = sorted([f'src/{f}' for f in os.listdir('src') if cohort in f and pheno in f and '.regions.manhattan' in f])
            qq_files = sorted([f'src/{f}' for f in os.listdir('src') if cohort in f and pheno in f and '.regions.qq' in f])

            for manhattan_file, qq_file in zip(manhattan_files, qq_files):
                img(src=manhattan_file, cls='img-fluid')
                br()
                img(src=qq_file, cls='img-fluid')
                br()
            
            hits_file = 'src/saige_exwas_suggestive_regions.csv'
            hits = pd.read_csv(hits_file).rename(columns=region_col_map).rename(columns={'Pvalue': 'P'})
            hits = hits.sort_values(by='P', ascending=True)
            hits = hits[(hits['COHORT'] == cohort) & (hits['PHENO'] == pheno)]
            keepCols = ['Region', 'Group', 'max_MAF', 'BETA_Burden', 'SE_Burden', 'P']
            hits = hits[keepCols].rename(columns={'imputationInfo': 'IMPUTE_SCORE'}).drop_duplicates()
            makeHitsTable(hits)

    with div(id=singles_html_id, cls='collapse multi-collapse'):
        with div(cls='card card-body'):
            h3('SAIGE-GENE Singles:')
            manhattan_file = 'src/{c}.{p}.singles.manhattan_vertical.png'.format(c=cohort, p=pheno)
            img(src=manhattan_file, cls='img-fluid')
            br()
            qq_file = 'src/{c}.{p}.singles.qq.png'.format(c=cohort, p=pheno)
            img(src=qq_file, cls='img-fluid')
            br()
            hits_file = 'src/saige_exwas_suggestive_singles.csv'
            hits = pd.read_csv(hits_file).rename(columns=singles_col_map).rename(columns={'p.value': 'P'})
            hits = hits.sort_values(by='P', ascending=True)
            keepCols = ['MarkerID', 'Allele1', 'Allele2', 'AF_Allele2', 'BETA', 'SE', 'P']
            hits = hits[(hits['COHORT'] == cohort) & (hits['PHENO'] == pheno)]
            hits = hits[keepCols].rename(columns={'imputationInfo': 'IMPUTE_SCORE'}).drop_duplicates()
            makeHitsTable(hits)


def makeSingleQuantitativeReport(cohort, pheno):
    # Get Source Files (hits table and plots)
    src_files = ['src/' + f for f in os.listdir('src') if cohort in f and pheno in f]
    if len(src_files) < 2:
        return
    if counts.loc[(cohort, pheno), 'N'] < 500:
        return
    
    p('Mean: ' + str(counts.loc[(cohort, pheno), 'mean'].astype(int)), style='margin-left:10%')
    p('STD: ' + str(counts.loc[(cohort, pheno), 'std'].astype(int)), style='margin-left:10%')
    br()
    p('Min: ' + str(counts.loc[(cohort, pheno), 'min'].astype(int)), style='margin-left:10%')
    p('25th %ile: ' + str(counts.loc[(cohort, pheno), '25%'].astype(int)), style='margin-left:10%')
    p('Median: ' + str(counts.loc[(cohort, pheno), '50%'].astype(int)), style='margin-left:10%')
    p('75th %ile: ' + str(counts.loc[(cohort, pheno), '75%'].astype(int)), style='margin-left:10%')
    p('Max: ' + str(counts.loc[(cohort, pheno), 'max'].astype(int)), style='margin-left:10%')
    br()

    region_html_id = 'collapse' + '-'.join([cohort, 'regions'])
    singles_html_id = 'collapse' + '-'.join([cohort, 'singles'])

    with p():
        butt = button('SAIGE-GENE Region Results (' + cohort + ')', cls='btn btn-light')
        butt['type'] = 'button'
        butt['data-bs-toggle'] = 'collapse'
        butt['data-bs-target'] = '#' + region_html_id
        butt['aria-expanded'] = 'false'
        butt['aria-controls'] = region_html_id
        
        butt = button('SAIGE-GENE Singles Results (' + cohort + ')', cls='btn btn-light')
        butt['type'] = 'button'
        butt['data-bs-toggle'] = 'collapse'
        butt['data-bs-target'] = '#' + singles_html_id
        butt['aria-expanded'] = 'false'
        butt['aria-controls'] = singles_html_id

    with div(id=region_html_id, cls='collapse multi-collapse'):
        with div(cls='card card-body'):
            h3('SAIGE-GENE Regions:')
            manhattan_files = sorted([f'src/{f}' for f in os.listdir('src') if cohort in f and pheno in f and '.regions.manhattan' in f])
            qq_files = sorted([f'src/{f}' for f in os.listdir('src') if cohort in f and pheno in f and '.regions.qq' in f])

            for manhattan_file, qq_file in zip(manhattan_files, qq_files):
                img(src=manhattan_file, cls='img-fluid')
                br()
                img(src=qq_file, cls='img-fluid')
                br()
            
            hits_file = 'src/saige_exwas_suggestive_regions.csv'
            hits = pd.read_csv(hits_file).rename(columns=region_col_map).rename(columns={'Pvalue': 'P'})
            hits = hits.sort_values(by='P', ascending=True)
            hits = hits[(hits['COHORT'] == cohort) & (hits['PHENO'] == pheno)]
            keepCols = ['Region', 'Group', 'max_MAF', 'BETA_Burden', 'SE_Burden', 'P']
            hits = hits[keepCols].rename(columns={'imputationInfo': 'IMPUTE_SCORE'}).drop_duplicates()
            makeHitsTable(hits)

    with div(id=singles_html_id, cls='collapse multi-collapse'):
        with div(cls='card card-body'):
            h3('SAIGE-GENE Singles:')
            manhattan_file = 'src/{c}.{p}.singles.manhattan_vertical.png'.format(c=cohort, p=pheno)
            img(src=manhattan_file, cls='img-fluid')
            br()
            qq_file = 'src/{c}.{p}.singles.qq.png'.format(c=cohort, p=pheno)
            img(src=qq_file, cls='img-fluid')
            br()
            hits_file = 'src/saige_exwas_suggestive_singles.csv'
            hits = pd.read_csv(hits_file).rename(columns=singles_col_map).rename(columns={'p.value': 'P'})
            hits = hits.sort_values(by='P', ascending=True)
            hits = hits[(hits['COHORT'] == cohort) & (hits['PHENO'] == pheno)]
            keepCols = ['MarkerID', 'Allele1', 'Allele2', 'AF_Allele2', 'BETA', 'SE', 'P']
            hits = hits[keepCols].rename(columns={'imputationInfo': 'IMPUTE_SCORE'}).drop_duplicates()
            makeHitsTable(hits)

args = make_arg_parser().parse_args()
p_cutoff = args.pval
cohorts = args.cohorts
all_phenos = args.phenotypes
output_dir = args.outDir
region_colnames_file = args.regionColnames
region_colnames_rows = open(region_colnames_file).read().splitlines()
region_col_map = dict(zip([r.split('=')[1] for r in region_colnames_rows],
                   [r.split('=')[0] for r in region_colnames_rows]))


singles_colnames_file = args.singlesColnames
singles_colnames_rows = open(singles_colnames_file).read().splitlines()
singles_col_map = dict(zip([r.split('=')[1] for r in singles_colnames_rows],
                   [r.split('=')[0] for r in singles_colnames_rows]))


print(f'Phenotypes: {all_phenos}')
print(f'Cohorts: {cohorts}')
print(f'P-Value Cutoff: {p_cutoff}')

counts = pd.read_csv('src/pheno_summaries.csv', dtype={'PHENO':str})
counts = counts.set_index(['COHORT', 'PHENO'])

# Generate each sub-page
cohorts = sorted(counts.index.get_level_values('COHORT').unique())
if 'Prevalence' in counts.columns:
    bin_phenos = counts.dropna(subset='Prevalence').index.get_level_values('PHENO').unique()
else:
    bin_phenos = []

if 'mean' in counts.columns:
    quant_phenos = counts.dropna(subset='mean').index.get_level_values('PHENO').unique()
else:
    quant_phenos = []

# Generate index page
index_file = makeMainPage(bin_phenos, quant_phenos)
# specify output directory if given
if output_dir:
    outfile=f'{output_dir}/index.html'
else:
    outfile=f'index.html'
open(outfile, 'w+').write(str(index_file))

all_phenos = [p for p in bin_phenos]
all_phenos.extend(quant_phenos)

for pheno in all_phenos:
    codePage = dominate.document(title = pheno + ': ExWAS Results Report')
    with codePage.head:
        addCssJsToHead()
    with codePage:
        with div(cls='container'):
            header(a('Back to Index', href='index.html'))
            h1(pheno)
            if pheno in bin_phenos:
                img(src='src/{p}.barplots.png'.format(p=pheno), cls='img-fluid')
            else:
                img(src='src/{p}.violinplot.png'.format(p=pheno), cls='img-fluid')
            with div():
                with p():
                    for cohort in cohorts:
                        if pheno in bin_phenos and counts.loc[(cohort, pheno), 'Cases'] < 50:
                            continue
                        if pheno in quant_phenos and counts.loc[(cohort, pheno), 'N'] < 500:
                            continue
                        butt = button('Cohort: ' + cohort, cls='btn btn-light')
                        butt['type'] = 'button'
                        butt['data-bs-toggle'] = 'collapse'
                        butt['data-bs-target'] = '#' + 'collapse' + cohort
                        butt['aria-expanded'] = 'false'
                        butt['aria-controls'] = 'collapse' + cohort
                for cohort in cohorts:
                    with div(id='collapse' + cohort, cls='collapse multi-collapse'):
                        with div(cls='card card-body'):
                            if pheno in bin_phenos:
                                makeSingleBinaryReport(cohort, pheno)
                            else:
                                makeSingleQuantitativeReport(cohort, pheno)
    # specify output directory if given
    if output_dir:
        outfile=f'{output_dir}/{pheno}.html'
    else:
        outfile=f'{pheno}.html'
    open(outfile, 'w+').write(str(codePage))
