nextflow.enable.dsl = 2

include {SAIGE_GWAS} from './saige_gwas.nf'
include {GWAMA_META} from './gwama_meta.nf'

workflow {
    saige_output = SAIGE_GWAS()
    saige_sumstats = saige_output[0]
    meta_sumstats = GWAMA_META(saige_sumstats)
}