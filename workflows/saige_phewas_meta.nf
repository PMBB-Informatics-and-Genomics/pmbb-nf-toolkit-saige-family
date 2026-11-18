import groovy.json.JsonOutput
import org.yaml.snakeyaml.Yaml

params.biobank_yamls = []
params.analyses = [:]
params.allow_missing_cohorts = false
params.p_cutoff_summarize = 1e-5
params.min_meta_studies = 2
params.meta_p_col = 'p_ivw'
params.test_id_cols = ['exposure', 'outcome']
params.pheno_descriptions_file = null
params.forest_plot_max = 50
params.varphewas_plot_max_labels = 20
params.varphewas_plot_max_labels_per_category = 5
params.my_python = '/opt/conda/bin/python'
params.host = ""
params.skipPostprocessingErrors = false

include {
    harmonize_phewas_sumstats
    run_phewas_meta
    plot_phewas_meta
    forest_plots_phewas_meta
    make_phewas_meta_report
    testIdCols
} from '../processes/saige_phewas_meta.nf'

include {
    get_script_file_names
    dump_params_to_json
} from '../processes/saige_helpers.nf'

// ---------------------------------------------------------------------------------------------
// YAML resolution. A biobank YAML inherits one study YAML; see READMEs/PheWAS_Meta_docs.md.
// Every validation error stops the run before any task is submitted.
// ---------------------------------------------------------------------------------------------

String resolveAgainst(String path, String baseDir) {
    // Absolute paths are kept; relative ones resolve against baseDir
    def f = new File(path)
    return (f.isAbsolute() ? f : new File(baseDir, path)).toPath().normalize().toString()
}

Map loadYaml(String path) {
    def f = new File(path)
    if (!f.isFile()) {
        error "PheWAS-META: YAML file not found: ${path}"
    }
    def data
    try {
        data = new Yaml().load(f.text)
    } catch (Exception e) {
        error "PheWAS-META: could not parse YAML ${path}: ${e.message}"
    }
    if (!(data instanceof Map)) {
        error "PheWAS-META: ${path} must be a YAML mapping (key: value pairs)"
    }
    return data
}

Map stringMap(Object raw, String key, String path) {
    // YAML mapping -> LinkedHashMap<String, String> (YAML may parse names such as 008 or 250.2 as numbers)
    if (raw == null) {
        return [:]
    }
    if (!(raw instanceof Map)) {
        error "PheWAS-META: '${key}' in ${path} must be a mapping"
    }
    def out = new LinkedHashMap<String, String>()
    raw.each { k, v ->
        if (v == null) {
            error "PheWAS-META: '${key}' in ${path}: '${k}' has no value"
        }
        out[k.toString()] = v.toString()
    }
    return out
}

Map renameMap(Object raw, String key, String path) {
    def m = stringMap(raw, key, path)
    def dups = m.groupBy { it.value }.findAll { it.value.size() > 1 }
    if (dups) {
        def detail = dups.collect { canon, entries -> "${entries.keySet().join(', ')} -> ${canon}" }.join('; ')
        error "PheWAS-META: ${path}: '${key}' renames several local names to one canonical name: ${detail}"
    }
    return m
}

Map resolveBiobank(String yamlPath) {
    // Returns the biobank's merged config: its own keys plus exposures / outcomes / cont_exposure_transform
    // from the inherited study file. Relative inherits / results_dir resolve against the biobank YAML's folder.
    def path = new File(yamlPath).absoluteFile.toPath().normalize().toString()
    def dir = new File(path).parent
    def bb = loadYaml(path)

    def study_only = ['exposures', 'outcomes', 'cont_exposure_transform'].findAll { bb.containsKey(it) }
    if (study_only) {
        error "PheWAS-META: biobank YAML ${path} sets ${study_only.join(', ')}; these come only from the inherited study file (use 'rename' for local names)"
    }
    def allowed = ['inherits', 'biobank', 'results_dir', 'cohorts', 'rename']
    def unknown = bb.keySet().collect { it.toString() } - allowed
    if (unknown) {
        error "PheWAS-META: biobank YAML ${path} has unknown key(s) ${unknown.join(', ')} (allowed: ${allowed.join(', ')})"
    }
    def missing_keys = ['inherits', 'biobank', 'results_dir', 'cohorts'].findAll { bb[it] == null || bb[it].toString().trim() == '' }
    if (missing_keys) {
        error "PheWAS-META: biobank YAML ${path} is missing required key(s) ${missing_keys.join(', ')}"
    }
    def biobank = bb.biobank.toString()

    def study_path = resolveAgainst(bb.inherits.toString(), dir)
    def study = loadYaml(study_path)
    if (study.containsKey('inherits')) {
        error "PheWAS-META: study file ${study_path} (inherited by ${path}) has 'inherits'; a study file cannot inherit another file"
    }
    def study_allowed = ['exposures', 'outcomes', 'cont_exposure_transform']
    def study_unknown = study.keySet().collect { it.toString() } - study_allowed
    if (study_unknown) {
        error "PheWAS-META: study file ${study_path} has unknown key(s) ${study_unknown.join(', ')} (allowed: ${study_allowed.join(', ')})"
    }
    def study_missing = ['exposures', 'cont_exposure_transform'].findAll { study[it] == null || study[it].toString().trim() == '' }
    if (study_missing) {
        error "PheWAS-META: study file ${study_path} (inherited by ${path}) is missing required key(s) ${study_missing.join(', ')}"
    }
    def exposures = stringMap(study.exposures, 'exposures', study_path)
    if (!exposures) {
        error "PheWAS-META: study file ${study_path}: 'exposures' is empty"
    }
    def bad_types = exposures.findAll { !(it.value in ['binary', 'continuous']) }
    if (bad_types) {
        error "PheWAS-META: study file ${study_path}: exposure types must be 'binary' or 'continuous', got ${bad_types}"
    }
    def transform = study.cont_exposure_transform.toString()
    if (!(transform in ['none', 'inverse_normal'])) {
        error "PheWAS-META: study file ${study_path}: cont_exposure_transform must be 'none' or 'inverse_normal', got '${transform}'"
    }
    def outcomes = null
    if (study.outcomes != null) {
        if (!(study.outcomes instanceof List)) {
            error "PheWAS-META: study file ${study_path}: 'outcomes' must be a list"
        }
        outcomes = study.outcomes.collect { it.toString() }
        def dup_outcomes = outcomes.findAll { o -> outcomes.count(o) > 1 }.unique()
        if (dup_outcomes) {
            error "PheWAS-META: study file ${study_path}: duplicate outcome(s) ${dup_outcomes.join(', ')}"
        }
    }

    def cohorts = stringMap(bb.cohorts, 'cohorts', path)
    if (!cohorts) {
        error "PheWAS-META: biobank YAML ${path}: 'cohorts' is empty"
    }
    def rename = bb.rename == null ? [:] : bb.rename
    if (!(rename instanceof Map)) {
        error "PheWAS-META: biobank YAML ${path}: 'rename' must be a mapping with 'exposures' and/or 'outcomes'"
    }
    def rename_unknown = rename.keySet().collect { it.toString() } - ['exposures', 'outcomes']
    if (rename_unknown) {
        error "PheWAS-META: biobank YAML ${path}: unknown 'rename' key(s) ${rename_unknown.join(', ')} (allowed: exposures, outcomes)"
    }
    def rename_exposures = renameMap(rename.exposures, 'rename.exposures', path)
    def rename_outcomes = renameMap(rename.outcomes, 'rename.outcomes', path)
    rename_exposures.findAll { !exposures.containsKey(it.value) }.each { local, canon ->
        log.warn "PheWAS-META: ${path}: rename.exposures ${local} -> ${canon}: '${canon}' is not a study exposure, so this rename is unused"
    }

    // canonical exposure -> the name this biobank's run used (files and exposure_encoding are keyed by it)
    def inverse = rename_exposures.collectEntries { local, canon -> [(canon): local] }
    def local_exposures = exposures.keySet().collectEntries { canon -> [(canon): inverse.get(canon, canon)] }
    local_exposures.each { canon, local ->
        if (local == canon && rename_exposures.containsKey(local) && rename_exposures[local] != canon) {
            error "PheWAS-META: ${path}: local exposure '${local}' is renamed to '${rename_exposures[local]}', so it cannot also supply study exposure '${canon}'"
        }
    }

    def results_dir = resolveAgainst(bb.results_dir.toString(), dir)
    if (!new File(results_dir).isDirectory()) {
        error "PheWAS-META: biobank ${biobank} (${path}): results_dir not found: ${results_dir}"
    }
    def encoding = new File(results_dir, 'Exposure_Genotypes/exposure_encoding.tsv')
    if (!encoding.isFile()) {
        error "PheWAS-META: biobank ${biobank} (${path}): exposure encoding not found: ${encoding.path} (results_dir must be the launchDir of an SVP exposure-mode run)"
    }

    return [
        yaml: path, study_file: study_path, biobank: biobank, results_dir: results_dir, encoding: encoding.path,
        cohorts: cohorts, rename_exposures: rename_exposures, rename_outcomes: rename_outcomes,
        exposures: exposures, outcomes: outcomes, cont_exposure_transform: transform, local_exposures: local_exposures,
    ]
}

void checkStudiesAgree(List<Map> biobanks) {
    // Every biobank must inherit the same study content (and biobank names must be unique)
    def names = biobanks.collect { it.biobank }
    def dup_names = names.findAll { n -> names.count(n) > 1 }.unique()
    if (dup_names) {
        def detail = dup_names.collect { n -> "${n} (${biobanks.findAll { it.biobank == n }.collect { it.yaml }.join(', ')})" }.join('; ')
        error "PheWAS-META: biobank name(s) used by more than one YAML: ${detail}"
    }
    def ref = biobanks[0]
    biobanks.drop(1).each { b ->
        def diffs = []
        if (b.exposures != ref.exposures) {
            diffs << "exposures ${b.exposures} vs ${ref.exposures}"
        }
        if ((b.outcomes == null ? null : b.outcomes as Set) != (ref.outcomes == null ? null : ref.outcomes as Set)) {
            diffs << "outcomes ${b.outcomes} vs ${ref.outcomes}"
        }
        if (b.cont_exposure_transform != ref.cont_exposure_transform) {
            diffs << "cont_exposure_transform '${b.cont_exposure_transform}' vs '${ref.cont_exposure_transform}'"
        }
        if (diffs) {
            error "PheWAS-META: biobank ${b.biobank} (${b.yaml}, study ${b.study_file}) does not inherit the same study as ${ref.biobank} (${ref.yaml}, study ${ref.study_file}): ${diffs.join('; ')}"
        }
    }
}

void checkAnalysisLabels(Map analyses, List<Map> biobanks) {
    // Every analysis lists meta cohort labels defined exactly once across the biobank YAMLs
    if (!analyses) {
        error 'PheWAS-META: params.analyses is empty; map each analysis name to a list of meta cohort labels'
    }
    def name_ok = /[A-Za-z0-9_.-]+/
    def owners = [:]
    biobanks.each { b ->
        b.cohorts.each { run, label ->
            owners.get(label, []) << "${b.biobank} cohort ${run} (${b.yaml})"
        }
    }
    def bad_labels = owners.keySet().findAll { !(it ==~ name_ok) }
    if (bad_labels) {
        error "PheWAS-META: meta cohort label(s) ${bad_labels.join(', ')} may only contain letters, numbers, '_', '.' and '-'"
    }
    def dup_labels = owners.findAll { it.value.size() > 1 }
    if (dup_labels) {
        def detail = dup_labels.collect { label, where -> "${label} (${where.join(' and ')})" }.join('; ')
        error "PheWAS-META: meta cohort label(s) defined more than once: ${detail}"
    }
    analyses.each { name, labels ->
        if (!(name.toString() ==~ name_ok)) {
            error "PheWAS-META: analysis name '${name}' may only contain letters, numbers, '_', '.' and '-'"
        }
        if (!(labels instanceof List) || !labels) {
            error "PheWAS-META: analysis '${name}' must be a non-empty list of meta cohort labels, got ${labels}"
        }
        def ls = labels.collect { it.toString() }
        def undefined = ls.findAll { !owners.containsKey(it) }
        if (undefined) {
            error "PheWAS-META: analysis '${name}' uses undefined meta cohort label(s) ${undefined.join(', ')}; labels defined in the biobank YAMLs: ${owners.keySet().sort().join(', ')}"
        }
        def repeated = ls.findAll { l -> ls.count(l) > 1 }.unique()
        if (repeated) {
            error "PheWAS-META: analysis '${name}' lists meta cohort label(s) more than once: ${repeated.join(', ')}"
        }
    }
}

List harmonizeInputs(List<Map> biobanks, Set used_labels, boolean allow_missing) {
    // One [meta_cohort, config_json, existing sumstats, encoding] per biobank x run cohort used by an analysis.
    // Paths are built from config values (names may contain '__'); missing files go into the config's 'missing'.
    // A used run cohort with no sumstats file at all (usually a typo'd cohort name), or a study exposure with no
    // file in any used cohort, stops the run unless allow_missing (params.allow_missing_cohorts) is set.
    def inputs = []
    def found_exposures = [] as Set
    biobanks.each { b ->
        b.cohorts.each { run, label ->
            if (!(label in used_labels)) {
                log.info "PheWAS-META: ${b.biobank} cohort ${run} (${label}) is not in any analysis; skipped"
                return
            }
            def files = new LinkedHashMap<String, String>()
            def missing = new LinkedHashMap<String, String>()
            def staged = []
            b.local_exposures.each { canon, local ->
                def base = "${run}__${local}.chr1.variant_phewas.saige.gz".toString()
                def f = new File("${b.results_dir}/${run}__${local}/Sumstats/${base}")
                if (f.isFile()) {
                    files[local] = base
                    staged << file(f.path)
                    found_exposures << canon
                } else {
                    missing[local] = f.path
                    log.warn "PheWAS-META: ${b.biobank} cohort ${run}: no sumstats for exposure ${canon} (expected ${f.path}); skipped"
                }
            }
            if (!files) {
                def msg = "PheWAS-META: biobank ${b.biobank} run cohort '${run}' (meta cohort ${label}) matched no sumstats files in results_dir ${b.results_dir} (expected ${run}__<exposure>/Sumstats/${run}__<exposure>.chr1.variant_phewas.saige.gz); check the run cohort name in ${b.yaml}"
                if (!allow_missing) {
                    error "${msg}, or set allow_missing_cohorts = true to continue without it"
                }
                log.warn "${msg}; continuing because allow_missing_cohorts = true"
            }
            def config = [
                biobank: b.biobank, run_cohort: run, meta_cohort: label,
                cont_exposure_transform: b.cont_exposure_transform, exposures: b.exposures, outcomes: b.outcomes,
                rename_exposures: b.rename_exposures, rename_outcomes: b.rename_outcomes,
                files: files, missing: missing,
            ]
            inputs << [label, JsonOutput.toJson(config), staged, file(b.encoding)]
        }
    }
    def absent = biobanks[0].exposures.keySet().findAll { !(it in found_exposures) }
    if (absent) {
        def msg = "PheWAS-META: study exposure(s) ${absent.join(', ')} have no sumstats file in any used cohort of any biobank"
        if (!allow_missing) {
            error "${msg}; check the study and rename.exposures names, or set allow_missing_cohorts = true to continue without them"
        }
        log.warn "${msg}; continuing because allow_missing_cohorts = true"
    }
    return inputs
}

List yamlList(Object raw) {
    // params.biobank_yamls: a list, or (e.g. from --biobank_yamls a.yaml,b.yaml) a comma-separated string
    def items = raw == null ? [] : (raw instanceof List ? raw : [raw])
    return items.collectMany { it == null ? [] : it.toString().split(',').collect { s -> s.trim() } }.findAll { it }
}

workflow {
    // Resolve and validate the YAMLs before logging (and before any task), so a bad config fails fast
    def summary_dir = launchDir.toString() + '/Summary'
    def yamls = yamlList(params.biobank_yamls)
    if (!yamls) {
        error 'PheWAS-META: params.biobank_yamls is empty; list one YAML per biobank'
    }
    def biobanks = yamls.collect { resolveBiobank(resolveAgainst(it.toString(), launchDir.toString())) }
    checkStudiesAgree(biobanks)
    if (!(params.analyses instanceof Map)) {
        error "PheWAS-META: params.analyses must be a map of analysis name -> [meta cohort labels], got ${params.analyses}"
    }
    checkAnalysisLabels(params.analyses, biobanks)
    if (!(params.meta_p_col in ['p_ivw', 'p_stouffer'])) {
        error "PheWAS-META: meta_p_col must be 'p_ivw' or 'p_stouffer', got '${params.meta_p_col}'"
    }
    if (testIdCols() != ['exposure', 'outcome']) {
        error "PheWAS-META: test_id_cols must be ['exposure', 'outcome'], got ${params.test_id_cols}; variant mode (other test id columns) is not yet supported"
    }
    def used_labels = params.analyses.values().flatten().collect { it.toString() } as Set
    def allow_missing = params.allow_missing_cohorts.toString().toLowerCase() == 'true'
    // checks that each used run cohort has sumstats files (logs a warning per missing file)
    def harmonize_inputs = harmonizeInputs(biobanks, used_labels, allow_missing)

    log.info([
        "  NEXTFLOW - DSL2 - SAIGE PheWAS-META - P I P E L I N E",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "run as", workflow.commandLine),
        String.format("  %-25s : %s", "run location", launchDir),
        String.format("  %-25s : %s", "started at", workflow.start),
        String.format("  %-25s : %s", "my_python", params.my_python),
        String.format("  %-25s : %s", "host", params.host),
        String.format("  %-25s : %s", "skipPostprocessingErrors", params.skipPostprocessingErrors),
        "",
        "  Inputs and Analyses",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "biobank_yamls", params.biobank_yamls),
        String.format("  %-25s : %s", "analyses", params.analyses),
        String.format("  %-25s : %s", "allow_missing_cohorts", params.allow_missing_cohorts),
        String.format("  %-25s : %s", "study file", biobanks[0].study_file),
        String.format("  %-25s : %s", "study exposures", biobanks[0].exposures),
        String.format("  %-25s : %s", "study outcomes", biobanks[0].outcomes ?: 'all present'),
        "",
        "  Meta-Analysis Parameters",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "p_cutoff_summarize", params.p_cutoff_summarize),
        String.format("  %-25s : %s", "min_meta_studies", params.min_meta_studies),
        String.format("  %-25s : %s", "meta_p_col", params.meta_p_col),
        String.format("  %-25s : %s", "test_id_cols", params.test_id_cols),
        "",
        "  Plot Parameters",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "pheno_descriptions_file", params.pheno_descriptions_file),
        String.format("  %-25s : %s", "forest_plot_max", params.forest_plot_max),
        String.format("  %-25s : %s", "varphewas_plot_max_labels", params.varphewas_plot_max_labels),
        String.format("  %-25s : %s", "varphewas_plot_max_labels_per_category", params.varphewas_plot_max_labels_per_category),
        ""
    ].join("\n"))
    biobanks.each { b ->
        log.info "  ${b.biobank}: ${b.results_dir} cohorts ${b.cohorts} renames ${b.rename_exposures + b.rename_outcomes}"
    }

    def script_name_dict = get_script_file_names()

    // meta cohort -> biobank / run cohort for the report (only cohorts used by an analysis)
    def cohort_rows = biobanks.collectMany { b ->
        b.cohorts.findAll { run, label -> label in used_labels }.collect { run, label -> "${label},${b.biobank},${run}\n".toString() }
    }.sort()
    def cohorts_csv = Channel.fromList(cohort_rows).collectFile(
        name: 'phewas_meta_cohorts.csv', seed: 'meta_cohort,biobank,run_cohort\n', sort: false, storeDir: summary_dir)

    harmonize_phewas_sumstats(
        Channel.fromList(harmonize_inputs),
        file(script_name_dict['phewas_meta_harmonize'])
    )

    // [analysis, cohort_order, harmonized files in cohort_order]: each analysis keeps only its own labels
    def analyses = params.analyses.collect { name, labels -> [name.toString(), labels.collect { it.toString() }] }
    def meta_inputs = harmonize_phewas_sumstats.out.harmonized
        .toList()
        .flatMap { pairs ->
            def by_label = pairs.collectEntries { label, f -> [(label): f] }
            analyses.collect { name, labels -> [name, labels, labels.collect { by_label[it] }] }
        }
    run_phewas_meta(meta_inputs, file(script_name_dict['phewas_meta_tests']))

    // PheWAS plots (need outcome descriptions) and forest plots of the top hits
    if (params.pheno_descriptions_file) {
        plot_phewas_meta(
            run_phewas_meta.out.plot_input,
            file(script_name_dict['phewas_meta_phewas_plots']),
            file(params.pheno_descriptions_file)
        )
    }
    def cohort_orders = analyses.collectEntries { name, labels -> [(name): labels] }
    def forest_in = run_phewas_meta.out.filtered_by_analysis
        .join(run_phewas_meta.out.effects)
        .map { analysis, filtered, effects -> [analysis, cohort_orders[analysis], filtered, effects] }
    forest_plots_phewas_meta(forest_in, file(script_name_dict['phewas_meta_forest_plots']))

    // Summaries: each per-task CSV has a fixed header (written even when it has no rows)
    def suggestive_csv = run_phewas_meta.out.filtered.collectFile(
        name: 'phewas_meta_all_suggestive.csv', keepHeader: true, skip: 1, sort: { it.name }, storeDir: summary_dir)
    def test_counts_csv = run_phewas_meta.out.test_counts.collectFile(
        name: 'phewas_meta_test_counts.csv', keepHeader: true, skip: 1, sort: { it.name }, storeDir: summary_dir)
    def sample_sizes_csv = run_phewas_meta.out.sample_sizes.collectFile(
        name: 'phewas_meta_sample_sizes.csv', keepHeader: true, skip: 1, sort: { it.name }, storeDir: summary_dir)
    def harmonization_log_csv = harmonize_phewas_sumstats.out.log.collectFile(
        name: 'phewas_meta_harmonization_log.csv', keepHeader: true, skip: 1, sort: { it.name }, storeDir: summary_dir)

    dump_params_to_json(params, 'saige_phewas_meta')

    // HTML report (needs the PheWAS plots, so only with pheno_descriptions_file, as in Single-Variant PheWAS)
    if (params.pheno_descriptions_file) {
        def all_plots = plot_phewas_meta.out[1].mix(forest_plots_phewas_meta.out).flatten().collect().ifEmpty([])
        make_phewas_meta_report(
            all_plots,
            suggestive_csv,
            test_counts_csv,
            sample_sizes_csv,
            harmonization_log_csv,
            cohorts_csv,
            dump_params_to_json.out,
            file(script_name_dict['phewas_meta_manifest']),
            file(script_name_dict['phewas_meta_report'])
        )
    }
}
