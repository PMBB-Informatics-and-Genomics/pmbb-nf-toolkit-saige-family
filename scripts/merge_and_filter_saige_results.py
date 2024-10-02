import pandas as pd
import argparse as ap
import sys
import re

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")
    
    parser.add_argument('-s', '--sumstats', nargs='+', help='List of results files')
    parser.add_argument('-c', '--colnames', required=True, help='File with column name mappings')
    parser.add_argument('--cohort', help='Cohort', required=True)
    parser.add_argument('--pvalue', help='P-Value for top hits table', type=float)
    parser.add_argument('--casecontrol',help='Case Control Threshold',type=int,required=False)

    parser.add_argument('--chr', help='Chromosome')
    parser.add_argument('-p', '--pheno', help='Phenotype')

    parser.add_argument('--phewas', action='store_true', default=False, help='Use this flag if merging PheWAS results')
    parser.add_argument('--regions', action='store_true', default=False, help='Use this flag if merging gene group results')
    parser.add_argument('--rareVars', action='store_true', default=False, help='Use this flag if merging single rare variant ExWAS results')
    parser.add_argument('--counts', required=False, help='Phenotype summaries table so we can fill in sample size')
    
    return parser

def make_output_file_names(is_phewas, is_gene_burden, is_rare_vars, arg_dict):
    
    output_filenames = {}

    if not is_phewas and is_gene_burden and not is_rare_vars:
        # ExWAS
        pheno = arg_dict.pheno
        cohort = arg_dict.cohort
        output_filenames['merge'] = f'{cohort}.{pheno}.exwas_regions.saige.gz'
        output_filenames['filter'] = f'{cohort}.{pheno}.exwas_regions.filtered.saige.csv'
        output_filenames['merge_cauchy'] = f'{cohort}.{pheno}.exwas_cauchy_tests.saige.gz'
        output_filenames['merge_filter'] = f'{cohort}.{pheno}.exwas_cauchy_tests.filtered.saige.csv'
    elif not is_phewas and not is_gene_burden and not is_rare_vars:
        # GWAS
        pheno = arg_dict.pheno
        cohort = arg_dict.cohort

        output_filenames['merge'] = f'{cohort}.{pheno}.gwas.saige.gz'
        output_filenames['filter'] = f'{cohort}.{pheno}.gwas.filtered.saige.csv'
    elif is_phewas and is_gene_burden and not is_rare_vars:
        # Gene burden PheWAS
        chrom = arg_dict.chr
        cohort = arg_dict.cohort

        output_filenames['merge'] = f'{cohort}.chr{chrom}.gene_phewas_regions.saige.gz'
        output_filenames['filter'] = f'{cohort}.chr{chrom}.gene_phewas_regions.filtered.saige.csv'
        output_filenames['merge_cauchy'] = f'{cohort}.chr{chrom}.gene_phewas_cauchy_tests.saige.gz'
        output_filenames['merge_filter'] = f'{cohort}.chr{chrom}.gene_phewas_cauchy_tests.filtered.saige.csv'
    elif is_phewas and not is_gene_burden and not is_rare_vars:
        # Variant PheWAS
        chrom = arg_dict.chr
        cohort = arg_dict.cohort
        output_filenames['merge'] = f'{cohort}.chr{chrom}.variant_phewas.saige.gz'
        output_filenames['filter'] = f'{cohort}.chr{chrom}.variant_phewas.filtered.saige.csv'
    elif not is_phewas and not is_gene_burden and is_rare_vars:
        # Single rare-variant tests accompanying ExWAS
        pheno = arg_dict.pheno
        cohort = arg_dict.cohort
        output_filenames['merge'] = f'{cohort}.{pheno}.exwas_singles.saige.gz'
        output_filenames['filter'] = f'{cohort}.{pheno}.exwas_singles.filtered.saige.csv'
        output_filenames['merge_UR'] = f'{cohort}.{pheno}.exwas_ultrarare_tests.saige.gz'
        output_filenames['filter_UR'] = f'{cohort}.{pheno}.exwas_ultrarare_tests.filtered.saige.csv'
    elif is_phewas and not is_gene_burden and is_rare_vars:
        # Single rare-variant tests accompanying Gene burden PheWAS
        chrom = args.chr
        cohort = arg_dict.cohort
        output_filenames['merge'] = f'{cohort}.chr{chrom}.gene_phewas_singles.saige.gz'
        output_filenames['filter'] = f'{cohort}.chr{chrom}.gene_phewas_singles.filtered.saige.csv'
        output_filenames['merge_UR'] = f'{cohort}.chr{chrom}.gene_phewas_ultrarare_tests.saige.gz'
        output_filenames['filter_UR'] = f'{cohort}.chr{chrom}.gene_phewas_ultrarare_tests.filtered.saige.csv'
    
    return output_filenames

def make_column_name_map_dict(colnames_file):
    colnames_rows = open(colnames_file).read().splitlines()
    print(colnames_rows)
    col_map = dict(zip([r.split('=')[0] for r in colnames_rows],
                    [r.split('=')[1] for r in colnames_rows]))
    print(col_map)
    return(col_map)

ALL_SAIGE_COLS = ['phenotype', 'Region', 'Group', 'max_MAF', 'Pvalue', 
                  'Pvalue_Burden', 'Pvalue_SKAT', 'BETA_Burden', 
                  'SE_Burden', 'MAC', 'MAC_case', 'MAC_control', 
                  'Number_rare', 'Number_ultra_rare',
                  'CHR', 'POS', 'MarkerID', 'Allele1', 'Allele2', 'AC_Allele2',
                  'AF_Allele2', 'MissingRate', 'BETA', 'SE', 'Tstat', 'var', 'p.value',
                  'p.value.NA', 'Is.SPA', 'AF_case', 'AF_ctrl', 'N_case', 'N_ctrl',
                  'N_case_hom', 'N_case_het', 'N_ctrl_hom', 'N_ctrl_het', 'N']

def make_output_groups(is_gene_burden, is_rare_vars, all_results, output_file_dict):
    """
    Returns list of tuples with (DataFrame, unfiltered file name, filtered file name)
    """
    output_groups = []

    if is_gene_burden and not is_rare_vars:
        # ExWAS and Gene burden PheWAS results
        not_cauchy = all_results[all_results['Group'] != 'Cauchy']
        output_groups.append((not_cauchy, output_file_dict['merge'], output_file_dict['filter']))

        is_cauchy = all_results[all_results['Group'] == 'Cauchy']
        output_groups.append((is_cauchy, output_file_dict['merge_cauchy'], output_file_dict['merge_filter']))

    elif not is_gene_burden and not is_rare_vars:
        # GWAS and Variant PheWAS results
        output_groups.append((all_results, output_file_dict['merge'], output_file_dict['filter']))

    elif not is_gene_burden and is_rare_vars:
        # Single rare-variant tests accompanying ExWAS and Gene burden PheWAS
        not_UR = all_results[all_results['CHR'] != 'UR']
        output_groups.append((not_UR, output_file_dict['merge'], output_file_dict['filter']))

        is_UR = all_results[all_results['CHR'] == 'UR']
        output_groups.append((is_UR, output_file_dict['merge_UR'], output_file_dict['filter_UR']))
    
    return output_groups

args = make_arg_parser().parse_args()
colnames_file = args.colnames
input_files = args.sumstats
cohort = args.cohort
p_thresh = args.pvalue
cc_thresh = args.casecontrol

is_phewas = args.phewas
is_gene_burden = args.regions
is_rare_vars = args.rareVars

output_file_dict = make_output_file_names(is_phewas, is_gene_burden, is_rare_vars, args)
col_map = make_column_name_map_dict(args.colnames)

dfs = []
for f in input_files:
    print(f)
    # parse the file name to get the name of the phenotype
    if args.pheno is None:
        f_no_gz = f.replace('.gz', '')
        f_no_suffix = f_no_gz.replace('.singleAssoc', '').replace('.txt', '')
        f_no_chr = '.'.join(f_no_suffix.split('.')[:-1])
        f_no_cohort = f_no_chr.replace(cohort, '')
        f_pheno = f_no_cohort.split('/')[-1]
    else:
        f_pheno = args.pheno

    temp = pd.read_table(f, compression='gzip') # Read results file with default behaviors
    temp = temp.drop_duplicates() # The singleAssoc output repeats variants sometimes
    temp['phenotype'] = f_pheno # Add phenotype name to individual DataFrame
    dfs.append(temp) # Append this DataFrame to dfs

# Concatenate all DataFrames together
all = pd.concat(dfs)

if is_gene_burden:
    counts_table = pd.read_csv(args.counts)
    counts_table = counts_table[counts_table['COHORT'] == cohort].set_index('PHENO')
    all[['N_case', 'N_ctrl', 'N']] = counts_table.reindex(all['phenotype'])[['Cases', 'Controls', 'N']].values

# Make sure all results have total N
if 'N' not in all.columns:
    all['N'] = all[['N_case', 'N_ctrl']].sum(axis=1)

# Figure out the set of SAIGE columns that are in the results files
use_cols = [c for c in ALL_SAIGE_COLS if c in all.columns]
all = all[use_cols]
print(all)
print(all.columns)

# Get list of tuples with (DataFrame, Unfiltered File Name, Filtered File Name)
output_groups = make_output_groups(is_gene_burden, is_rare_vars, all, output_file_dict)

# Iterate over output groups
for df, out, out_filtered in output_groups:
    # Rename columns
    df = df.rename(columns=col_map)

    # Write unfiltered output
    df.to_csv(out, sep='\t', index=False, na_rep='NA')

    # Prepare filtered output
    p_col = 'Pvalue' if is_gene_burden else 'p.value'
    mapped_p_col = col_map[p_col] if p_col in col_map.keys() else p_col
    p_value_col = pd.to_numeric(df[mapped_p_col], errors='coerce')

    if cc_thresh:
    
        case_col = 'MAC_case' if is_gene_burden else 'N_case'
        mapped_case_col = col_map[case_col] if case_col in col_map.keys() else case_col
        case_col_value = pd.to_numeric(df[mapped_case_col], errors='coerce')

        control_col = 'MAC_control' if is_gene_burden else 'N_ctrl'
        mapped_control_col = col_map[control_col] if control_col in col_map.keys() else control_col
        control_col_value = pd.to_numeric(df[mapped_control_col], errors='coerce')
        df_filtered = df[(p_value_col <= p_thresh) & (case_col_value >= cc_thresh) & (control_col_value >= cc_thresh) ].copy()
        
        
    else:

        df_filtered = df[(p_value_col <= p_thresh)].copy()
        
    df_filtered.insert(0, 'cohort', cohort)
    # Write filtered output
    df_filtered.to_csv(out_filtered, index=False, na_rep='NA')