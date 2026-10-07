import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse

# read in args provided by nextflow channel to process calling python script
def make_arg_parser():
    parser = argparse.ArgumentParser(description="inputs for making single variant plots")
    
    parser.add_argument("--pheno_descriptions",help="file containinging phenotype name category and description")
    parser.add_argument("--input_file_list",nargs='+', help='List of results files')
    parsed = parser.parse_args()
    
    return parsed


# Function to read and concatenate files
def concatenate_files(files):
    dfs = []
    for file in files:
        df = pd.read_csv(file, compression='gzip', sep='\t')  # Assuming files are gzipped and tab-separated
        dfs.append(df)
    concatenated_df = pd.concat(dfs, ignore_index=True)
    return concatenated_df


args = make_arg_parser()
input_files =  args.input_file_list
pheno_descriptions_file_path = args.pheno_descriptions

# turn string rep of list into list
print(input_files)
files = [x.strip(',[]') for x in input_files]
# Concatenate files
concatenated_df = concatenate_files(files)


concatenated_df['p_value'] = concatenated_df['p_value'].astype('float')
filtered_df = concatenated_df[concatenated_df['p_value']<= .01]
pheno_descriptions_file = pd.read_csv(pheno_descriptions_file_path)
merged_df =pd.merge(concatenated_df, pheno_descriptions_file, left_on='phenotype', right_on='PHENO',how="left")
merged_df = merged_df.drop("DESCRIPTION",axis=1)
filtered_df = merged_df[merged_df['p_value']<= .01]


def read_phewas_results(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df

def plot_manhattan_with_annotations_for_variant(phewas_df, variant_id, phenotype_column, p_value_column, significance_threshold=0.05, save_path=None):

    variant_df = phewas_df[phewas_df['variant_id'] == variant_id]

    unique_phenotypes = variant_df[phenotype_column].unique()
    color_palette = sns.color_palette("husl", n_colors=len(unique_phenotypes))  # Use seaborn color palette
    color_map = {phenotype: color_palette[i] for i, phenotype in enumerate(unique_phenotypes)}

    fig, ax = plt.subplots(figsize=(12, 6))
    phenotype_indices = np.arange(len(unique_phenotypes))
    for i, phenotype in enumerate(unique_phenotypes):
        phenotype_df = variant_df[variant_df[phenotype_column] == phenotype]
        neg_log_p_values = -np.log10(phenotype_df[p_value_column])
        ax.scatter(phenotype_indices[i] + phenotype_df.index / len(variant_df), neg_log_p_values, color=color_map[phenotype], label=phenotype, alpha=0.7, s=10)


    ax.axhline(-np.log10(significance_threshold), color='gray', linestyle='--', linewidth=1)


    ax.set_xticks(phenotype_indices)
    ax.set_xticklabels(unique_phenotypes, rotation=45, ha='right')
    ax.set_xlabel('Phenotype')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'Manhattan Plot for Variant {variant_id}')
    ax.legend(title='Phenotype', bbox_to_anchor=(1.05, 1), loc='upper left')


    ax.set_xlim(-0.5, len(unique_phenotypes) - 0.5)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    

variant_ids = merged_df['variant_id'].unique()
for variant_id in variant_ids:
  filename = "Plotting.log.txt"
  plot_manhattan_with_annotations_for_variant(merged_df, variant_id, 'CATEGORY', 'p_value', significance_threshold=0.01, save_path=f'manhattan_plot_variant_{variant_id}.png')
  with open(filename, 'a') as file:  # Open the file in append mode
        file.write(f"Plot generated for variant_{variant_id} at manhattan_plot_variant_{variant_id}.png")
