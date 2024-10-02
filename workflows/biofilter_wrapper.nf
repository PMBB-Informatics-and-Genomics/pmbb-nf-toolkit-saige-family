
workflow {
    test_pos_file = "${launchDir}/test_positions_input.txt"
    test_pos_channel = Channel.of(new Tuple('test_pos', test_pos_file))
    test_output = BIOFILTER_POSITIONS(test_pos_channel)
}

workflow BIOFILTER_POSITIONS {
    take:
        data_positions // channel with tuples of (data nickname, position file)
    main:

        bf_annot = call_biofilter_positions(
            data_positions,
            params['biofilter_script'],
            params['biofilter_loki']
            )

        output_annot = assign_positions_rsids_genes(bf_annot)
    emit:
        output_annot
}

ANNOTATIONS = 'position_label snp position gene upstream downstream'

process call_biofilter_positions {
    publishDir "${launchDir}/Annotations/"

    input:
        tuple val(data_nickname), path(positions_file)
        path(biofilter_script)
        path(biofilter_loki)
    output:
        tuple val(data_nickname), path("${data_nickname}_biofilter_positions_annotations.txt")
    script:
        output_prefix = "${data_nickname}_biofilter_annotations"
        output_ext = ANNOTATIONS.split()[0] + '.' + ANNOTATIONS.split()[1..-1].join('-')
        """
        ${params.my_python} ${biofilter_script} \
          --verbose \
          --knowledge ${biofilter_loki} \
          --position-file ${positions_file} \
          --annotate ${ANNOTATIONS} \
          --report-invalid-input \
          --overwrite \
          --prefix ${output_prefix} \
          --ucsc-build-version ${params['biofilter_build']}

        mv ${output_prefix}.${output_ext} ${data_nickname}_biofilter_positions_annotations.txt
        """
    stub:
        """
        touch ${data_nickname}_biofilter_positions_annotations.txt
        """
}

process assign_positions_rsids_genes {
    publishDir "${launchDir}/Annotations"

    input:
        tuple val(data_nickname), path(annot_file)
    output:
        tuple val(data_nickname), path("${data_nickname}_biofilter_genes_rsids.csv")
    script:
        """
        #! ${params.my_python}

        import pandas as pd

        def assign_gene_annotations(biofilter_df, close_dist=5E4):
            # Picking first RSID by default because in dbSNP it's usually an SNV
            biofilter_df = biofilter_df[~biofilter_df.index.duplicated(keep='first')]

            # Define "close" as 50kb
            biofilter_df[['u-close', 'd-close']] = biofilter_df[['distance', 'distance.1']] < float(close_dist)
            biofilter_df['gene-combo'] = biofilter_df[['upstream', 'downstream']] \
                                        .apply(lambda x: '/'.join(x.astype(str)),axis=1)
            biofilter_df['missing'] = pd.isnull(biofilter_df['gene'])
            biofilter_df['use-combo'] = (biofilter_df['u-close'] == biofilter_df['d-close']) \
                                         & biofilter_df['missing']  # Both close or both far
            biofilter_df['use-u'] = ~biofilter_df['use-combo'] & biofilter_df['missing'] & biofilter_df['u-close']
            biofilter_df['use-d'] = ~biofilter_df['use-combo'] & biofilter_df['missing'] & biofilter_df['d-close']

            biofilter_df['Gene'] = biofilter_df['gene']
            for boolCol, geneCol in {'use-u': 'upstream', 'use-d': 'downstream', 'use-combo': 'gene-combo'}.items():
                idx = biofilter_df.index[biofilter_df[boolCol]]
                biofilter_df.loc[idx, 'Gene'] = biofilter_df.loc[idx, geneCol]

            return biofilter_df

        df = pd.read_table('${annot_file}', index_col='position_label')

        df = assign_gene_annotations(df, ${params['biofilter_close_dist']})

        outDF = df.reset_index()[['position_label', 'snp', 'chr', 'pos', 'Gene']]
        outDF = outDF.rename(columns={'position_label': 'Var_ID', 'snp': 'RSID',
                                      'chr': '#CHROM', 'pos': 'POS'})

        outDF.to_csv('${data_nickname}_biofilter_genes_rsids.csv', index=False)
        """
    stub:
        """
        touch ${data_nickname}_biofilter_genes_rsids.csv
        """
}
