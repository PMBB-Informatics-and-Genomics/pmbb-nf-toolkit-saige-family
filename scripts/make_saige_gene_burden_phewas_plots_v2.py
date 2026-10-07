# Imports
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse as ap

from adjustText import adjust_text

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    # Add a non-optional list argument for cohort
    parser.add_argument('-c', '--cohort', required=True, help='cohort to sample')
    parser.add_argument('-r', '--chromosome', required=True, help="Chromosome from process")
    parser.add_argument('-s', '--sumstats_file', required=True)
    parser.add_argument('-d', '--descriptions_file', required=True, help="Phenotype descriptions File")
    parser.add_argument('-p', '--pvalue_column', default=None, help="P-value column to use for plotting. Default: Regex ^(pval|p-value|p_value) on sumstats file")
    parser.add_argument('-g', '--gene_region_column', default=None, help="Gene Region Column used for subsetting. Default: Regex (gene|region|variant) on sumstats file")
    parser.add_argument('-b', '--beta_column', default=None, help="Beta column for plotting. Default: Regex ^(beta) on sumstats file")
    parser.add_argument('-m', '--maf_column', default=None, help="MAF column for subsetting unique values. Default: Regex ^(maf|max_maf) on sumstats file")
    parser.add_argument('-a', '--annotation_column', default=None, help="Group Annotation column for subsetting. Default: Regex ^(group|annot) on sumstats file")
    parser.add_argument('-l', '--gene_label', default=None, help="Label Column for Gene Region. Default: None")
    parser.add_argument('-x', '--pheno_col_sumstats', default=None, help="Phenotype ID Column for sumstats file. Default: Regex ^(phecode|pheno) on sumstats columns")
    parser.add_argument('-y', '--pheno_col_descriptions', default=None, help="Phenotype ID Column for descriptions file. Default: Regex ^(phecode|pheno) on descriptions columns")

    return parser

def read_table(file_path):
    # Determine the delimiter based on the file extension
    # compression = "infer"
    if file_path.endswith('.csv'):
        delimiter = ','
    elif file_path.endswith('.tsv'):
        delimiter = '\t'
    elif file_path.endswith('.gz'):
        delimiter = '\t'
    else:
        raise ValueError("Unsupported file extension. Please provide a CSV or TSV file.")
    
    # Read the file using pandas
    df = pd.read_csv(file_path, delimiter=delimiter)
    return df

def get_category_counts(df):
    """
    Function to get category counts for subdf

    Returns:
        _type_: _description_
    """
    counts = df.groupby('category').size()
    counts_dict = counts.to_dict()

    # sort dfframe by counts
    categories_order = counts.sort_values(ascending=False).index.tolist()
    df.loc[:,'category'] = pd.Categorical(df['category'].copy(deep=True), categories=categories_order)
    df = df.sort_values(by = "category")

    # make labels
    label_map = {}
    for key in counts_dict:
        label_map[key] = f"{key} ({counts_dict[key]})"
    return label_map

def check_and_prepend(value):
    '''
    Function to check and modify the phecode values if needed
    '''
    if not value.lower().startswith('phe'):
        return 'Phe' + value
    else:
        return value
    
def sort_df_categories(df):
    """
    Function to sort df by category counts

    Returns:
        _type_: _description_
    """
    counts = df.groupby('category').size()

    # sort dfframe by counts
    categories_order = counts.sort_values(ascending=False).index.tolist()
    df.loc[:,'category'] = pd.Categorical(df['category'], categories=categories_order)
    df = df.sort_values(by = "category")
    return df
    
def select_top_percent(group):
    """
    Function to select the top 10% of values within each group

    Args:
        group (_type_): _description_

    Returns:
        _type_: _description_
    """
    n_top = max(int(len(group) * 0.1), 1)  # Calculate the top 10%, ensuring at least one
    return group.nsmallest(n_top, 'LOG10P')

# Parse script arguments
args = make_arg_parser().parse_args()
sumstats_file = args.sumstats_file
cohort = args.cohort
gene_region_column = args.gene_region_column
descriptions_file = args.descriptions_file
pvalue_column = args.pvalue_column
maf_column = args.maf_column
beta_column = args.beta_column
annotation_column = args.annotation_column
pheno_col_sumstats = args.pheno_col_sumstats
pheno_col_descriptions = args.pheno_col_descriptions
chromosome = args.chromosome


# Read in summary stats file
sumstats = read_table(sumstats_file)

# read in descriptions file
pheno_descriptions = pd.read_csv(descriptions_file)

# find phecode columns in the data files
phecode_pattern = r'^(phecode|pheno)'
if not pheno_col_sumstats:
    pheno_col_sumstats = sumstats.columns[sumstats.columns.str.match(phecode_pattern, case=False)].tolist()[0]

if not pheno_col_descriptions:
    pheno_col_descriptions = pheno_descriptions.columns[pheno_descriptions.columns.str.match(phecode_pattern, case=False)].tolist()[0]

# Apply the function to the phecode column to prepend Phe if needed
# sumstats[pheno_col_sumstats] = sumstats[pheno_col_sumstats].astype(str).apply(check_and_prepend)
# pheno_descriptions[pheno_col_descriptions] = pheno_descriptions[pheno_col_descriptions].astype(str).apply(check_and_prepend)

# get descriptions and category columns
descriptions_pattern = r'^(descr)'
description_col = pheno_descriptions.columns[pheno_descriptions.columns.str.match(descriptions_pattern, case=False)].tolist()[0]
category_pattern = r'^(cat)'
category_col = pheno_descriptions.columns[pheno_descriptions.columns.str.match(category_pattern, case=False)].tolist()[0]
pheno_descriptions.rename(columns={category_col: "category"},inplace=True)
pheno_descriptions.rename(columns={description_col: "description"},inplace=True)
# fill in unknown and make own category
pheno_descriptions['category'] = pheno_descriptions['category'].fillna("Unknown")
pheno_descriptions['description'] = pheno_descriptions['description'].fillna(pheno_descriptions['category'])
# merge sumstats and descriptions
data = pd.merge(left=sumstats, right=pheno_descriptions, left_on=pheno_col_sumstats, right_on=pheno_col_descriptions)

# Identify column names to be used later
# --------------------------------------
# find gene region column
if not gene_region_column:
    region_pattern = r'^(gene|region)$'
    gene_region_column = data.columns[data.columns.str.match(region_pattern, case=False)].tolist()[0]

# find p-value column
if not pvalue_column:
    # if id doesn't exist, return first match to pvalue substring search
    pval_pattern = r'^(pval|p-value|p_value)'
    pvalue_column = data.columns[data.columns.str.match(pval_pattern, case=False)].tolist()[0]

# find MAF column
if not maf_column:
    maf_pattern = r'^(maf|max_maf)'
    maf_column = data.columns[data.columns.str.match(maf_pattern, case=False)].tolist()[0]
    
# find Beta column
if not beta_column:
    beta_pattern = r'^(beta)'
    beta_column = data.columns[data.columns.str.match(beta_pattern, case=False)].tolist()[0]

# find Group column
if not annotation_column:
    annotation_pattern = r'^(group|annot)'
    annotation_column = data.columns[data.columns.str.match(annotation_pattern, case=False)].tolist()[0]

# print(f"Region Column: {gene_region_column}\nPVAL Column: {pvalue_column}\nMAF Column: {maf_column}\nGroup Column: {annotation_column}\n")

# remove nas in case something weird happened
check_cols = [gene_region_column,pvalue_column,maf_column, beta_column,annotation_column]
data.dropna(subset=check_cols,inplace=True)

# Log-transform the p values for plotting
data['LOG10P'] = -np.log10(data[pvalue_column])

# Loop over 
for region, subDF1 in data.groupby(gene_region_column): # Iterate over genes
    for group, subDF2 in subDF1.groupby(annotation_column): # Iterate over annotation groups
        for maf, subDF3 in subDF2.groupby(maf_column): # Iterate over max MAF thresholds
            output_file = f'{cohort}.{region}.{group.replace(";","-")}.maf{str(maf).replace(".","-")}.phewas.png'
            plot_title = f'PheWAS Manhattan for {region} in {cohort}\nMax MAF = {maf}, Annotation Group = {group}'
            print(f"\n{plot_title}\n-->{output_file}")
            # subset data based on region, cohort, and maf
            subdf = data[(data[gene_region_column] == region) & (data[annotation_column] == group) & (data[maf_column] == maf)]
            # subdf = sort_df_categories(subdf).reset_index(drop=True)
            # get label map with counts
            label_map = get_category_counts(subdf)
            subdf['category_label'] = subdf['category'].map(label_map)
            # Get and Plot Bonferroni corrected P Value Threshold
            p_thresh = -np.log10(0.05 / len(subdf))
            # Filter points with pvalues greater than 0.005 for each unique label
            filtered_points = subdf[subdf['LOG10P'] > p_thresh]
            # Find top hit (highest pvalue) for each label1 if no points with pvalue > 0.005
            top_hits = subdf.loc[subdf.groupby('category_label',observed=True)['LOG10P'].idxmax()]
            # Merge filtered points and top hits
            labeled_points = pd.concat([filtered_points, top_hits]).drop_duplicates()
            
            # Filter points retaining only top 10 percent for each unique label
            # labeled_points = subdf.groupby('category_label',group_keys=False).apply(select_top_percent)
            
            # make the plot
            # ------------------------------------------------------------------
            # Get the color palette
            palette = sns.color_palette('deep',subdf.category_label.nunique())
            
            # get order
            counts = subdf.groupby('category_label').size()
            label_order = counts.sort_values(ascending=False).index.tolist()

            # do the plot
            g = sns.catplot(subdf, x="category_label", y="LOG10P", height=10,
                            aspect=1.5, palette=palette, legend=False,
                            s=30, jitter = 0.5, order = label_order)

            # Label the points with the same color as the plot
            # create empty list to hold text objects
            palette_dict = dict(zip(label_order,palette))
            txt_plots = []
            for i, (_, row) in enumerate(labeled_points.iterrows()):
                # get the color of the datapoint
                color = palette_dict[row['category_label']]
                # create a plot_text object and append to list
                txt = plt.text(x = row['category_label'], y = row['LOG10P'], 
                s=row['description'], ha='center', va='bottom', color=color,size=10,
                bbox=dict(facecolor='white', edgecolor='black', linewidth=0,boxstyle='round,pad=.2',alpha=0.5)
                )
                txt_plots.append(txt)

            # use adjustTEXT package to nudge the text
            adjust_text(txt_plots)

            # add horizontal line for thresh (catplot) if there are points above
            if len(filtered_points) > 1:
                g.refline(y = p_thresh, color = "gray", lw = 1)

            # set title and axis
            g.set(title=plot_title, xlabel='PheCode Category (n)', ylabel='-log10p')

            # set rotation
            g.set_xticklabels(rotation=45)
            plt.tight_layout()
            
            # save
            plt.savefig(output_file)


sumstats_file

def plots_filename(row):
    cohort = row['cohort']
    region = row[gene_region_column]
    ann = row[annotation_column]
    maf = row[maf_column]
    output_plot_file = f'Plots/{cohort}.{region}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.phewas.png'
    # print(f"Debug output plots file: \n{output_plot_file}\n")
    return output_plot_file

data.insert(0,'cohort',cohort)
grouped_df = data.groupby(['cohort',gene_region_column, annotation_column,maf_column]).first().reset_index()
grouped_df['results_file'] = f"{cohort}/Sumstats/{sumstats_file}"
grouped_df['gene_phewas_plot'] = grouped_df.apply(lambda row: plots_filename(row), axis=1)
# grouped_df['gene_phewas_plot'] = grouped_df.apply(plots_filename, axis=1)
grouped_df['filtered_results_file'] = "Summary/saige_exwas_suggestive_regions.csv"
output_file = f"{cohort}.{chromosome}.phewas_regions.plots_manifest.csv"
grouped_df.to_csv(output_file,index=False)