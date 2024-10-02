import pandas as pd
import numpy as np
import dominate
from dominate.tags import *
from dominate.util import raw
import os
import sys
import argparse as ap

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")
    
    # Add a non-optional list argument for phenotypes
    parser.add_argument('-p', '--phenotypes', nargs='+', required=True, help='List of phenotypes')
    
    # Add a non-optional list argument for cohorts
    parser.add_argument('-c', '--cohorts', nargs='+', required=True, help='List of cohorts')

    # Add arguments for regions and singles column name files
    parser.add_argument('--regionsColnames')
    parser.add_argument('--singlesColnames')

    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default=None, help='Path to output directory. Default: current working directory')
    
    # Add an argument for p-value cutoff
    parser.add_argument('--pval', required=True, help='P-Value Cutoff for Hit Tables', type=float)
    
    return parser


def addCssJsToHead():
    '''
    Adds bootstrap CSS resources to the HTML header so we can use their styling elements
    '''
    link(href='https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/css/bootstrap.min.css',
        rel='stylesheet',
        integrity='sha384-4bw+/aepP/YC94hEpVNVgiZdgIC5+VKNBQNGCHeKRQN+PtmoHDEXuppvnDJzQIu9',
        crossorigin='anonymous')
    script(src='https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/js/bootstrap.bundle.min.js',
            integrity='sha384-HwwvtgBNo3bZJJLYd8oVXjrBZt8cqVSpeBNS5n7C8IVInixGAoxmnlMuBnhbgrkm',
            crossorigin='anonymous')
    meta(charset='utf-8', name='viewport', content='width=device-width')


def makeMainPage(bin_phenos, quant_phenos):
    '''
    Constructs the main page of the report with links to individual phenotype pages
    :param bin_phenos: list of binary phenotypes used in this analysis
    :type bin_phenos: list
    :param quant_phenos: list of binary phenotypes used in this analysis
    :type quant_phenos: list
    :return: main page dominate HTML document object
    :rtype: dominate.document()
    '''

    mainPage = dominate.document(title = 'ExWAS Results Reports')

    # Add CSS to header
    with mainPage.head:
        addCssJsToHead()
    
    # Add link to methods description
    with mainPage:
        a('Methods:', href='src/saige_exwas_methods.html')

    # Add binary phenotype links as needed
    if len(bin_phenos) > 0:
        with mainPage:
            h2('Binary Phenotypes:')
            for pheno in bin_phenos:
                a(pheno, href=pheno + '.html')
                br()

    # Add quantitative phenotype links as needed
    if len(quant_phenos) > 0:
        with mainPage:
            h2('Quantitative Phenotypes:')
            for pheno in quant_phenos:
                a(pheno, href=pheno + '.html')
                br()

    return mainPage


def makeHitsTable(hits):
    '''
    Creates an HTML table with the top hits from the analysis
    :param hits: a DataFrame with the results from one cohort and phenotype
    :type hits: pd.DataFrame()
    '''
    # If the DataFrame is empty, add a message
    if (len(hits) == 0):
        p('No Hits Met User-Specified P-Value Filter ({p:.2E})'.format(p=p_cutoff), style='color:red;margin-left:10%')
        return
    
    # Copy the DataFrame subset
    printHits = hits.copy()
    # Note the user-specified filter for reference
    p('User-Specified P-Value Filter: {p:.2E}'.format(p=p_cutoff), style='color:red;margin-left:10%')
    # Format the P-values with scientific notation
    printHits['P'] = printHits['P'].apply(lambda x: '{:.2E}'.format(x))

    # Nested blocks to create the table:
    # Div > table > tr > th or td
    with div(cls='table-responsive'):
        with table(cls="table table-striped"):
            # Create Table Header Row
            with tr():
                for c in hits.columns:
                    th(c, style='border: 2px solid black')
            
            # Create Table Data Rows
            for _, row in printHits.iterrows():
                with tr():
                    for _, item in row.items():
                        try:
                            td(str(np.round(item, 5)), style='border: 2px solid black')
                        except:
                            td(item, style='border: 2px solid black')


def makeSingleBinaryReport(cohort, pheno):
    '''
    Generates a part of a phenotype page for one cohort
    :param cohort: The population tested in this section of results
    :type cohort: str
    :param pheno: The phenotype tested in this section of results
    :type pheno: str
    '''

    # Get Source Files (hits table and plots)
    src_files = ['src/' + f for f in os.listdir('src') if cohort in f and pheno in f]
    if len(src_files) < 2:
        return
    if counts.loc[(cohort, pheno), 'Cases'] < 50:
        return
    
    # Start with some basic info about the cohort-phenotype combo
    p('Cases: ' + str(counts.loc[(cohort, pheno), 'Cases'].astype(int)), style='margin-left:10%')
    p('Controls: ' + str(counts.loc[(cohort, pheno), 'Controls'].astype(int)), style='margin-left:10%')

    # Generate HTML IDs to use with the collapsible elements
    region_html_id = 'collapse' + '-'.join([cohort, 'regions'])
    singles_html_id = 'collapse' + '-'.join([cohort, 'singles'])

    # Add buttons for the collapsible elements
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

    # Add Regions Results
    with div(id=region_html_id, cls='collapse multi-collapse'):
        with div(cls='card card-body'):
            h3('SAIGE-GENE Regions:')

            # Images
            manhattan_files = sorted([f'src/{f}' for f in os.listdir('src') if cohort in f and pheno in f and '.regions.manhattan' in f])
            qq_files = sorted([f'src/{f}' for f in os.listdir('src') if cohort in f and pheno in f and '.regions.qq' in f])

            manhattan_file_df = pd.DataFrame([f.replace('src/', '').split('.') for f in manhattan_files], index=manhattan_files)
            print(manhattan_file_df)
            manhattan_file_df.columns = ['COHORT', 'PHENO', 'ANNOT', 'MAF', 'REGIONS', 'PLOT', 'EXT']

            maf_filters = list(manhattan_file_df['MAF'].unique())
            annot_filters = list(manhattan_file_df['ANNOT'].unique())

            with form(id=f"optionsForm={cohort}-{pheno}"):
                p('MAF')
                for i, maf in enumerate(maf_filters):
                    with label(maf):
                        input_(type='radio', name=f'maf-{cohort}-{pheno}', value=maf, checked= i == 0)
                p('Annot')
                for i, annot in enumerate(annot_filters):
                    with label(annot):
                        input_(type='radio', name=f'annot-{cohort}-{pheno}', value=annot, checked= i== 0)
                
                img(id=f'selectedImage{cohort}{pheno}', src='', cls='img-fluid')
                img(id=f'selectedImage{cohort}{pheno}QQ', src='', cls='img-fluid')

                raw(f"""
                    <script>
                    // Define the base URL for the images
                    const baseURL{cohort}{pheno} = 'src/{cohort}.{pheno}';

                    // Get references to the radio buttons and the image element
                    const maf{cohort}{pheno}Radios = document.getElementsByName('maf-{cohort}-{pheno}');
                    const annot{cohort}{pheno}Radios = document.getElementsByName('annot-{cohort}-{pheno}');
                    const selectedImage{cohort}{pheno} = document.getElementById('selectedImage{cohort}{pheno}');

                    // Function to update the image based on the selected options
                    function updateImage{cohort}{pheno}() {{
                    // Get the selected MAF and ANNOT values
                    const selectedMaf = document.querySelector('input[name="maf-{cohort}-{pheno}"]:checked').value;
                    const selectedAnnot = document.querySelector('input[name="annot-{cohort}-{pheno}"]:checked').value;

                    // Construct the image URL
                    const imageURL = `${{baseURL{cohort}{pheno}}}.${{selectedAnnot}}.${{selectedMaf}}.regions.manhattan.png`;
                    console.log(imageURL)

                    // Update the image source
                    selectedImage{cohort}{pheno}.src = imageURL;
                    }}

                    // Function to update the image based on the selected options
                    function updateImage{cohort}{pheno}QQ() {{
                    // Get the selected MAF and ANNOT values
                    const selectedMaf = document.querySelector('input[name="maf-{cohort}-{pheno}"]:checked').value;
                    const selectedAnnot = document.querySelector('input[name="annot-{cohort}-{pheno}"]:checked').value;

                    // Construct the image URL
                    const imageURL = `${{baseURL{cohort}{pheno}}}.${{selectedAnnot}}.${{selectedMaf}}.regions.qq.png`;
                    console.log(imageURL)

                    // Update the image source
                    selectedImage{cohort}{pheno}QQ.src = imageURL;
                    }}

                    // Add event listeners to the radio buttons to trigger the image update
                    maf{cohort}{pheno}Radios.forEach(radio => radio.addEventListener('change', updateImage{cohort}{pheno}));
                    annot{cohort}{pheno}Radios.forEach(radio => radio.addEventListener('change', updateImage{cohort}{pheno}));
                    maf{cohort}{pheno}Radios.forEach(radio => radio.addEventListener('change', updateImage{cohort}{pheno}QQ));
                    annot{cohort}{pheno}Radios.forEach(radio => radio.addEventListener('change', updateImage{cohort}{pheno}QQ));

                    // Initially update the image based on the default selected options
                    updateImage{cohort}{pheno}();
                    updateImage{cohort}{pheno}QQ();
                    </script>
                """)

            """
            for manhattan_file, qq_file in zip(manhattan_files, qq_files):
                img(src=manhattan_file, cls='img-fluid')
                br()
                img(src=qq_file, cls='img-fluid')
                br()
            """
            
            # Hits Table
            hits_file = 'src/saige_exwas_suggestive_regions.csv'
            hits = pd.read_csv(hits_file).rename(columns=region_col_map).rename(columns={'Pvalue': 'P'})
            hits = hits.sort_values(by='P', ascending=True)
            hits = hits[(hits['cohort'] == cohort) & (hits['phenotype'] == pheno)]
            keepCols = ['Region', 'Group', 'max_MAF', 'BETA_Burden', 'SE_Burden', 'P']
            hits = hits[keepCols].rename(columns={'imputationInfo': 'IMPUTE_SCORE'}).drop_duplicates()
            makeHitsTable(hits)

    # Add Singles Results
    with div(id=singles_html_id, cls='collapse multi-collapse'):
        with div(cls='card card-body'):
            h3('SAIGE-GENE Singles:')

            # Images
            manhattan_file = 'src/{c}.{p}.singles.manhattan_vertical.png'.format(c=cohort, p=pheno)
            img(src=manhattan_file, cls='img-fluid')
            br()
            qq_file = 'src/{c}.{p}.singles.qq.png'.format(c=cohort, p=pheno)
            img(src=qq_file, cls='img-fluid')
            br()

            # Hits Table
            hits_file = 'src/saige_exwas_suggestive_singles.csv'
            hits = pd.read_csv(hits_file).rename(columns=singles_col_map).rename(columns={'p.value': 'P'})
            hits = hits.sort_values(by='P', ascending=True)
            keepCols = ['MarkerID', 'Allele1', 'Allele2', 'AF_Allele2', 'BETA', 'SE', 'P']
            hits = hits[(hits['cohort'] == cohort) & (hits['phenotype'] == pheno)]
            hits = hits[keepCols].rename(columns={'imputationInfo': 'IMPUTE_SCORE'}).drop_duplicates()
            makeHitsTable(hits)


def makeSingleQuantitativeReport(cohort, pheno):
    '''
    Generates a part of a phenotype page for one cohort
    :param cohort: The population tested in this section of results
    :type cohort: str
    :param pheno: The phenotype tested in this section of results
    :type pheno: str
    '''

    # Get Source Files (hits table and plots)
    src_files = ['src/' + f for f in os.listdir('src') if cohort in f and pheno in f]
    if len(src_files) < 2:
        return
    if counts.loc[(cohort, pheno), 'N'] < 500:
        return
    
    # Add some basic summary of the cohort-phenotype combo
    p('Mean: ' + str(counts.loc[(cohort, pheno), 'mean'].astype(int)), style='margin-left:10%')
    p('STD: ' + str(counts.loc[(cohort, pheno), 'std'].astype(int)), style='margin-left:10%')
    br()
    p('Min: ' + str(counts.loc[(cohort, pheno), 'min'].astype(int)), style='margin-left:10%')
    p('25th %ile: ' + str(counts.loc[(cohort, pheno), '25%'].astype(int)), style='margin-left:10%')
    p('Median: ' + str(counts.loc[(cohort, pheno), '50%'].astype(int)), style='margin-left:10%')
    p('75th %ile: ' + str(counts.loc[(cohort, pheno), '75%'].astype(int)), style='margin-left:10%')
    p('Max: ' + str(counts.loc[(cohort, pheno), 'max'].astype(int)), style='margin-left:10%')
    br()

    # HTML IDs for the collapsible elements
    region_html_id = 'collapse' + '-'.join([cohort, 'regions'])
    singles_html_id = 'collapse' + '-'.join([cohort, 'singles'])

    # Buttons for the collapsible elements
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

    # Add Regions Results
    with div(id=region_html_id, cls='collapse multi-collapse'):
        with div(cls='card card-body'):
            h3('SAIGE-GENE Regions:')

            # Images
            manhattan_files = sorted([f'src/{f}' for f in os.listdir('src') if cohort in f and pheno in f and '.regions.manhattan' in f])
            qq_files = sorted([f'src/{f}' for f in os.listdir('src') if cohort in f and pheno in f and '.regions.qq' in f])

            manhattan_file_df = pd.DataFrame([f.replace('src/', '').split('.') for f in manhattan_files], index=manhattan_files)
            print(manhattan_file_df)
            manhattan_file_df.columns = ['COHORT', 'PHENO', 'ANNOT', 'MAF', 'REGIONS', 'PLOT', 'EXT']

            maf_filters = list(manhattan_file_df['MAF'].unique())
            annot_filters = list(manhattan_file_df['ANNOT'].unique())

            with form(id=f"optionsForm={cohort}-{pheno}"):
                p('MAF')
                for i, maf in enumerate(maf_filters):
                    with label(maf):
                        input_(type='radio', name=f'maf-{cohort}-{pheno}', value=maf, checked= i == 0)
                p('Annot')
                for i, annot in enumerate(annot_filters):
                    with label(annot):
                        input_(type='radio', name=f'annot-{cohort}-{pheno}', value=annot, checked= i== 0)
                
                img(id=f'selectedImage{cohort}{pheno}', src='', cls='img-fluid')
                img(id=f'selectedImage{cohort}{pheno}QQ', src='', cls='img-fluid')

                raw(f"""
                    <script>
                    // Define the base URL for the images
                    const baseURL{cohort}{pheno} = 'src/{cohort}.{pheno}';

                    // Get references to the radio buttons and the image element
                    const maf{cohort}{pheno}Radios = document.getElementsByName('maf-{cohort}-{pheno}');
                    const annot{cohort}{pheno}Radios = document.getElementsByName('annot-{cohort}-{pheno}');
                    const selectedImage{cohort}{pheno} = document.getElementById('selectedImage{cohort}{pheno}');

                    // Function to update the image based on the selected options
                    function updateImage{cohort}{pheno}() {{
                    // Get the selected MAF and ANNOT values
                    const selectedMaf = document.querySelector('input[name="maf-{cohort}-{pheno}"]:checked').value;
                    const selectedAnnot = document.querySelector('input[name="annot-{cohort}-{pheno}"]:checked').value;

                    // Construct the image URL
                    const imageURL = `${{baseURL{cohort}{pheno}}}.${{selectedAnnot}}.${{selectedMaf}}.regions.manhattan.png`;
                    console.log(imageURL)

                    // Update the image source
                    selectedImage{cohort}{pheno}.src = imageURL;
                    }}

                    // Function to update the image based on the selected options
                    function updateImage{cohort}{pheno}QQ() {{
                    // Get the selected MAF and ANNOT values
                    const selectedMaf = document.querySelector('input[name="maf-{cohort}-{pheno}"]:checked').value;
                    const selectedAnnot = document.querySelector('input[name="annot-{cohort}-{pheno}"]:checked').value;

                    // Construct the image URL
                    const imageURL = `${{baseURL{cohort}{pheno}}}.${{selectedAnnot}}.${{selectedMaf}}.regions.qq.png`;
                    console.log(imageURL)

                    // Update the image source
                    selectedImage{cohort}{pheno}QQ.src = imageURL;
                    }}

                    // Add event listeners to the radio buttons to trigger the image update
                    maf{cohort}{pheno}Radios.forEach(radio => radio.addEventListener('change', updateImage{cohort}{pheno}));
                    annot{cohort}{pheno}Radios.forEach(radio => radio.addEventListener('change', updateImage{cohort}{pheno}));
                    maf{cohort}{pheno}Radios.forEach(radio => radio.addEventListener('change', updateImage{cohort}{pheno}QQ));
                    annot{cohort}{pheno}Radios.forEach(radio => radio.addEventListener('change', updateImage{cohort}{pheno}QQ));

                    // Initially update the image based on the default selected options
                    updateImage{cohort}{pheno}();
                    updateImage{cohort}{pheno}QQ();
                    </script>
                """)

            """
            for manhattan_file, qq_file in zip(manhattan_files, qq_files):
                img(src=manhattan_file, cls='img-fluid')
                br()
                img(src=qq_file, cls='img-fluid')
                br()
            """
            
            # Hits Table
            hits_file = 'src/saige_exwas_suggestive_regions.csv'
            hits = pd.read_csv(hits_file).rename(columns=region_col_map).rename(columns={'Pvalue': 'P'})
            hits = hits.sort_values(by='P', ascending=True)
            hits = hits[(hits['cohort'] == cohort) & (hits['phenotype'] == pheno)]
            keepCols = ['Region', 'Group', 'max_MAF', 'BETA_Burden', 'SE_Burden', 'P']
            hits = hits[keepCols].rename(columns={'imputationInfo': 'IMPUTE_SCORE'}).drop_duplicates()
            makeHitsTable(hits)

    # Add Singles Results
    with div(id=singles_html_id, cls='collapse multi-collapse'):
        with div(cls='card card-body'):
            h3('SAIGE-GENE Singles:')

            # Images
            manhattan_file = 'src/{c}.{p}.singles.manhattan_vertical.png'.format(c=cohort, p=pheno)
            img(src=manhattan_file, cls='img-fluid')
            br()
            qq_file = 'src/{c}.{p}.singles.qq.png'.format(c=cohort, p=pheno)
            img(src=qq_file, cls='img-fluid')
            br()

            # Hits Table
            hits_file = 'src/saige_exwas_suggestive_singles.csv'
            hits = pd.read_csv(hits_file).rename(columns=singles_col_map).rename(columns={'p.value': 'P'})
            hits = hits.sort_values(by='P', ascending=True)
            hits = hits[(hits['cohort'] == cohort) & (hits['phenotype'] == pheno)]
            keepCols = ['MarkerID', 'Allele1', 'Allele2', 'AF_Allele2', 'BETA', 'SE', 'P']
            hits = hits[keepCols].rename(columns={'imputationInfo': 'IMPUTE_SCORE'}).drop_duplicates()
            makeHitsTable(hits)

# Parse Script Arguments
args = make_arg_parser().parse_args()
p_cutoff = args.pval
cohorts = args.cohorts
all_phenos = args.phenotypes
output_dir = args.outDir

# Turn the column name files into dictionaries
region_colnames_file = args.regionsColnames
region_colnames_rows = open(region_colnames_file).read().splitlines()
region_col_map = dict(zip([r.split('=')[1] for r in region_colnames_rows],
                   [r.split('=')[0] for r in region_colnames_rows]))

singles_colnames_file = args.singlesColnames
singles_colnames_rows = open(singles_colnames_file).read().splitlines()
singles_col_map = dict(zip([r.split('=')[1] for r in singles_colnames_rows],
                   [r.split('=')[0] for r in singles_colnames_rows]))


# Helpful script-level information
print(f'Phenotypes: {all_phenos}')
print(f'Cohorts: {cohorts}')
print(f'P-Value Cutoff: {p_cutoff}')

# Read phenotype counts
counts = pd.read_csv('src/pheno_summaries.csv', dtype={'phenotype':str})
counts = counts.set_index(['COHORT', 'PHENO'])

# Get lists of binary and quantitative phenotypes
cohorts = sorted(counts.index.get_level_values('COHORT').unique())
if 'Prevalence' in counts.columns:
    bin_phenos = counts.dropna(subset='Prevalence').index.get_level_values('PHENO').unique()
else:
    bin_phenos = []

if 'mean' in counts.columns:
    quant_phenos = counts.dropna(subset='mean').index.get_level_values('PHENO').unique()
else:
    quant_phenos = []

# Generate main index page
index_file = makeMainPage(bin_phenos, quant_phenos)

# specify output directory if given
if output_dir:
    outfile=f'{output_dir}/index.html'
else:
    outfile=f'index.html'
open(outfile, 'w+').write(str(index_file))

# Make list of all phenotypes
all_phenos = [p for p in bin_phenos]
all_phenos.extend(quant_phenos)

# Iterate over phenotypes to make individual phenotype report pages
for pheno in all_phenos:
    # Start new document
    codePage = dominate.document(title = pheno + ': ExWAS Results Report')
    # Add CSS to Header
    with codePage.head:
        addCssJsToHead()
    
    # Nested HTML elements start here:
    with codePage:
        with div(cls='container'):
            # General phenotype-level things for the top of the page
            # Div > header, h1, images for phenotype summary counts/distributions
            header(a('Back to Index', href='index.html'))
            h1(pheno)
            if pheno in bin_phenos:
                img(src='src/{p}.barplots.png'.format(p=pheno), cls='img-fluid')
            else:
                img(src='src/{p}.violinplot.png'.format(p=pheno), cls='img-fluid')
            
            # Cohort-specific sections
            # div > p > buttons
            # div > div > div (collapsible cards for reports)
            with div():
                with p():
                    for cohort in cohorts:
                        # Add collapsible buttons
                        butt = button('Cohort: ' + cohort, cls='btn btn-light')
                        butt['type'] = 'button'
                        butt['data-bs-toggle'] = 'collapse'
                        butt['data-bs-target'] = '#' + 'collapse' + cohort
                        butt['aria-expanded'] = 'false'
                        butt['aria-controls'] = 'collapse' + cohort
                
                # Add individual cohort reports
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
