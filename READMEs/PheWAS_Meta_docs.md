Documentation for PheWAS-META
============================

# Module Overview

PheWAS-META meta-analyzes **phenotype exposure PheWAS** results across cohorts and biobanks. Its inputs are the output folders of Single-Variant PheWAS runs in exposure mode (see [Phenotype Exposure Mode](Single-Variant_PheWAS_docs.md#phenotype-exposure-mode)). Each biobank runs SVP on its own data and shares the results folder. PheWAS-META then:

1. reads one YAML file per biobank, each inheriting a shared study YAML that fixes the canonical exposure and outcome names;
2. harmonizes every biobank × run cohort (renames local names, checks exposure types and the continuous-exposure transform, drops self-tests and outcomes outside the study list);
3. meta-analyzes each user-defined **analysis** (a list of meta cohorts) per exposure × outcome: fixed-effect inverse-variance (IVW), signed sample-size weighted Stouffer, and Cochran's Q / I²;
4. writes per-analysis tables and run-wide summary CSVs.

No genotype data or SAIGE runs are involved. Every step is a small Python task in the `pennbiobank/saige` image.

## Software Requirements

* [Nextflow version 24.04.3](https://www.nextflow.io/docs/latest/cli.html)

  > **Warning:** Nextflow 25.x and 26.x introduce breaking changes incompatible with this pipeline. Use version 24.04.3 exactly.

* [Singularity 3.8.3](https://sylabs.io/docs/) OR [Docker 4.30.0](https://docs.docker.com/)

## Commands for Running the Workflow

* Docker: `docker pull pennbiobank/saige:latest`; Singularity: `singularity build saige.sif docker://pennbiobank/saige:latest`
* Run from your analysis directory (the launch directory; all outputs are written there):

```sh
cd /path/to/analysis/phewas_meta_run1
nextflow run /path/to/toolkit/workflows/saige_phewas_meta.nf \
    -c nextflow.config \
    -c saige_phewas_meta.config \
    -profile standard
```

* Add `-resume` to pick up an interrupted run. `conf/saige_phewas_meta.config` is an example params file to copy and edit.

# Configuration

A run has three kinds of configuration:

| What | Where | Holds |
|------|-------|-------|
| Study YAML (one, shared) | e.g. `study.yaml` | canonical exposure names and types, optional outcome list, `cont_exposure_transform` |
| Biobank YAML (one per biobank) | e.g. `bb_a.yaml` | which study file it inherits, the SVP results folder, run cohort → meta cohort labels, local-name renames |
| Nextflow params | e.g. `saige_phewas_meta.config` | the list of biobank YAMLs, the analyses, thresholds and plot settings |

## Study YAML

```yaml
cont_exposure_transform: none            # every biobank's SVP run must have used this ('none' or 'inverse_normal')
exposures:                                # canonical name -> type, same words as exposure_encoding.tsv
  LDL_median: continuous
  smoker: binary
outcomes: [T2D, CAD, "250.2"]             # optional; omitted = every outcome present in the sumstats
```

| Key | Required | Description |
|-----|----------|-------------|
| `exposures` | **Required** | Mapping of canonical exposure name → `binary` or `continuous`. Each harmonize task checks this against the run's `Exposure_Genotypes/exposure_encoding.tsv` and fails on a mismatch. |
| `cont_exposure_transform` | **Required** | `none` or `inverse_normal`. Every continuous exposure in every run must have been encoded with this transform (the encoding's `transform` column); binary exposures must say `none`. A mismatch fails the run, because betas on different scales cannot be pooled. |
| `outcomes` | Optional | List of canonical outcome names to keep. Rows for other outcomes are dropped (and logged as `dropped_not_in_study`). |

A study file cannot itself have `inherits`, and no other keys are allowed.

## Biobank YAML

```yaml
inherits: study.yaml                      # relative to this file's folder
biobank: BBA
results_dir: /path/to/bba_svp_run         # launchDir of the biobank's SVP exposure-mode run; relative = this file's folder
cohorts:                                  # SVP run cohort -> meta cohort label (used in params.analyses)
  AFR_ALL: BBA_AFR
  EUR_ALL: BBA_EUR
rename:                                   # optional: local name -> canonical name
  exposures: {ldl_mgdl: LDL_median}
  outcomes: {t2d_local: T2D}
```

| Key | Required | Description |
|-----|----------|-------------|
| `inherits` | **Required** | Path to the study YAML. Relative paths resolve against the folder of this biobank YAML. |
| `biobank` | **Required** | Biobank name. Must be unique across `biobank_yamls`. Written to the `biobank` column of the harmonized tables and log. |
| `results_dir` | **Required** | The launch directory of the biobank's SVP exposure-mode run. Relative paths resolve against the folder of this biobank YAML. It must contain `Exposure_Genotypes/exposure_encoding.tsv`. |
| `cohorts` | **Required** | Mapping of SVP run cohort (an entry of that run's `cohort_list`) → meta cohort label. Labels must be unique across all biobanks and may contain only letters, numbers, `_`, `.` and `-`. |
| `rename.exposures` | Optional | Local exposure name (as used in that biobank's run) → canonical study name. Two local names may not map to one canonical name. |
| `rename.outcomes` | Optional | Local outcome name → canonical study name. Same rule. |

**Merge rule:** `exposures`, `outcomes` and `cont_exposure_transform` come **only** from the inherited study file. A biobank YAML that sets any of them is an error, and every biobank must inherit the same study content (the same exposures and types, the same outcome set and the same transform; the study files may live in different folders). Local naming differences go through the biobank's `rename` maps only.

**Quote names that look like numbers** (phecodes like `"008"` or `"250.10"`): YAML otherwise reads them as numbers, and `250.10` would become `250.1`.

### Which files are read

For every run cohort `C` in `cohorts` and every study exposure `E` (with local name `L`, i.e. `E` after reversing `rename.exposures`), the workflow expects:

```
<results_dir>/C__L/Sumstats/C__L.chr1.variant_phewas.saige.gz
<results_dir>/Exposure_Genotypes/exposure_encoding.tsv
```

Paths are built from the configured names and never split on `__`, so names may contain `__`. A missing sumstats file is logged (Nextflow warning plus a `missing_file` row in the harmonization log) and the run continues. But a used run cohort with **no** sumstats file for any exposure (usually a typo'd cohort name in `cohorts`), or a study exposure that has no file in any used cohort of any biobank, stops the run before any task starts, naming the biobank, run cohort and `results_dir`; set `allow_missing_cohorts = true` to log a warning and continue instead. A missing `exposure_encoding.tsv` or `results_dir` stops the run before any task starts. Run cohorts whose label is not used in any analysis are skipped.

**Sharing results from a remote biobank** (e.g. All of Us, where the SVP run cannot be read directly): only these files are needed. Copy them out keeping the same directory layout under a folder, and point that biobank's `results_dir` at the folder:

```
<results_dir>/Exposure_Genotypes/exposure_encoding.tsv
<results_dir>/<cohort>__<exposure>/Sumstats/<cohort>__<exposure>.chr1.variant_phewas.saige.gz   # one per run cohort x exposure
```

### Validation

These errors stop the run before any task is submitted, naming the offending file, key or label:

- a biobank YAML that sets `exposures`, `outcomes` or `cont_exposure_transform`, has an unknown key, or is missing `inherits`, `biobank`, `results_dir` or `cohorts`
- a study YAML that has `inherits` or an unknown key, is missing `exposures` or `cont_exposure_transform`, or uses an exposure type other than `binary`/`continuous`
- two local names renamed to one canonical name, or a local name that is renamed but is also needed as a study exposure under its own name
- biobanks that inherit different study content, or two YAMLs with the same `biobank` name
- a meta cohort label defined twice, an analysis that names an undefined label (or lists one twice), or an empty `analyses`
- a used run cohort with no sumstats files, or a study exposure with no sumstats in any used cohort (unless `allow_missing_cohorts = true`)
- `meta_p_col` other than `p_ivw`/`p_stouffer`, or `test_id_cols` other than `['exposure', 'outcome']`

These fail inside the harmonize task (the task's `ERROR:` message is printed in the Nextflow log):

- an exposure whose type or transform in `exposure_encoding.tsv` disagrees with the study
- a sumstats file whose `exposure` column is not the expected exposure
- duplicate (exposure, outcome) rows after renaming

## Nextflow parameters

```groovy
params {
    my_python = '/opt/conda/bin/python'
    biobank_yamls = ['/path/to/bb_a.yaml', '/path/to/bb_b.yaml']
    analyses = [
        'ALL': ['BBA_AFR', 'BBA_EUR', 'BBB_AFR', 'BBB_EUR'],
        'AFR': ['BBA_AFR', 'BBB_AFR'],
    ]
    p_cutoff_summarize = 1e-5
    min_meta_studies = 2
}
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `biobank_yamls` | list of paths | `[]` (**required**) | One biobank YAML per biobank. Relative paths resolve against the launch directory. On the command line, give a comma-separated list: `--biobank_yamls bb_a.yaml,bb_b.yaml`. |
| `analyses` | map: name → list of labels | `[:]` (**required**) | Each analysis meta-analyzes the listed meta cohort labels; the list order sets the order of the `direction` string and of the per-cohort columns in the effects table. Analysis names may contain only letters, numbers, `_`, `.` and `-`. |
| `allow_missing_cohorts` | bool | `false` | When `false`, a used run cohort with no sumstats file for any exposure (e.g. a typo'd cohort name), or a study exposure with no file in any used cohort, stops the run before any task. When `true`, these only log a warning (and `missing_file` rows) and the run continues without them. |
| `p_cutoff_summarize` | float | `1e-5` | A test is a top hit (written to `.filtered.csv` and `phewas_meta_all_suggestive.csv`) when `meta_p_col` ≤ this value... |
| `min_meta_studies` | int | `2` | ...and at least this many cohorts contributed (`N_studies`). Applies to top hits only; the full table keeps every test. |
| `meta_p_col` | string | `'p_ivw'` | P-value column used to sort and filter: `p_ivw` or `p_stouffer`. |
| `test_id_cols` | list of strings | `['exposure', 'outcome']` | Columns that define one meta-analysis test. Must be `['exposure', 'outcome']`; any other value stops the run (variant mode is not yet supported). On the command line, give it comma-separated (`--test_id_cols exposure,outcome`). |
| `pheno_descriptions_file` | path | `null` | CSV with columns `PHENO`, `DESCRIPTION`, `CATEGORY` used to label and group outcomes in PheWAS plots (same format as SVP). |
| `forest_plot_max` | int | `50` | Maximum number of forest plots (top hits by `meta_p_col`) per analysis. |
| `varphewas_plot_max_labels` | int | `20` | Maximum number of outcome labels on each PheWAS plot (same as SVP). |
| `varphewas_plot_max_labels_per_category` | int | `5` | Maximum number of labels per phenotype category on each PheWAS plot (same as SVP). |
| `my_python` | path | `'/opt/conda/bin/python'` | Python executable inside the container. |
| `host` | string | `""` | Execution environment hint. Set to `"AOU"` for All of Us (larger base memory for every step), or leave empty. |
| `skipPostprocessingErrors` | bool | `false` | When `true`, a failed plot step (PheWAS or forest plots) is ignored instead of stopping the run. |

Memory is requested dynamically: harmonize, meta-analysis and PheWAS plot tasks start at 8 GB (16 GB with `host = "AOU"`), forest plots at 4 GB (8 GB), and a task killed for memory (exit 137–140) is retried up to 3 times with the base × attempt. Exception: the PheWAS and forest plot tasks carry the `safe_to_skip` label, so when `conf/nextflow.config` is loaded its `safe_to_skip` error strategy replaces this retry — a plot task killed for memory is then skipped (or stops the run, per `skipPostprocessingErrors`) instead of retried; raise its memory with a `withName` selector if needed. Any other failure, such as a harmonization error, stops the run without a retry. For scale: meta-analysis of 300 exposures × 2,000 outcomes × 10 cohorts (6 M input rows) took about 100 s and 2 GB.

# Methods

For each analysis and each exposure × outcome, let cohort *i* have effect β<sub>i</sub>, standard error SE<sub>i</sub>, p-value p<sub>i</sub> and effective sample size n<sub>eff,i</sub>. Cohorts without that test are left out (`?` in the `direction` string).

- **Effective sample size:** `n_eff = 2 / (1/n_case + 1/n_ctrl)` for binary outcomes, `n` for quantitative outcomes.
- **IVW (fixed effect):** w<sub>i</sub> = 1/SE<sub>i</sub>², β<sub>IVW</sub> = Σw<sub>i</sub>β<sub>i</sub> / Σw<sub>i</sub>, SE<sub>IVW</sub> = √(1/Σw<sub>i</sub>), two-sided p from β<sub>IVW</sub>/SE<sub>IVW</sub>.
- **Signed Stouffer:** Z<sub>i</sub> = sign(β<sub>i</sub>) · Φ⁻¹(1 − p<sub>i</sub>/2), Z = Σ√n<sub>eff,i</sub> Z<sub>i</sub> / √(Σn<sub>eff,i</sub>), two-sided p. The sign matters: an unsigned Stouffer (as in ExWAS-META) adds effects in opposite directions together instead of letting them cancel. p-values of 0 are clipped at 1e-300.
- **Heterogeneity:** Cochran's Q = Σw<sub>i</sub>(β<sub>i</sub> − β<sub>IVW</sub>)² with inverse-variance weights, p from χ² with k−1 degrees of freedom, and I² = max(0, (Q − (k−1))/Q). Both are empty when only one cohort contributes.
- **Transform check:** continuous-exposure betas are per unit of the encoded exposure, so pooling is only valid if every run used the same `cont_exposure_transform`. The harmonize step enforces this from each run's `exposure_encoding.tsv`.
- **Units:** the workflow assumes a continuous exposure has the **same units** in every biobank (e.g. LDL in mg/dL everywhere). It does not convert units; harmonize them before the SVP runs.
- **Self-tests:** rows where the outcome equals the exposure (after renaming) are dropped.

# Outputs

All paths are relative to the launch directory.

| Path | Contents |
|------|----------|
| `Harmonized/<meta_cohort>.harmonized.tsv.gz` | Long table for one biobank × run cohort: `biobank, cohort, exposure, exposure_type, outcome, outcome_type, beta, standard_error, p_value, n, n_case, n_ctrl, n_eff` (canonical names; `cohort` is the meta cohort label). |
| `Harmonized/<meta_cohort>.harmonization_log.csv` | `biobank, run_cohort, meta_cohort, item_type, name, action, detail`. Actions: `renamed`, `missing_file`, `empty_file`, `dropped_self_test`, `dropped_not_in_study`, `dropped_invalid_stats`. |
| `Meta_Sumstats/<analysis>/<analysis>.phewas_meta.tsv.gz` | One row per exposure × outcome: `exposure_type, outcome_type, N_studies, direction, n_total, n_case_total, n_ctrl_total, n_eff_total, beta_ivw, se_ivw, p_ivw, z_stouffer, p_stouffer, q_het, p_het, i2_het`. Sorted by exposure, then `meta_p_col`. |
| `Meta_Sumstats/<analysis>/<analysis>.phewas_meta.filtered.csv` | Top hits (see `p_cutoff_summarize` and `min_meta_studies`), with an `analysis` column. May be header-only. |
| `Meta_Sumstats/<analysis>/<analysis>.phewas_meta.effects.tsv.gz` | Wide table: per-cohort `beta_`, `se_`, `p_`, `n_eff_<label>` columns in `analyses` order. |
| `Meta_Sumstats/<analysis>/<analysis>.phewas_meta.test_counts.csv` | Number of tests per exposure. |
| `Meta_Sumstats/<analysis>/<analysis>.phewas_meta.sample_sizes.csv` | Maximum `n_eff_total` and `N_studies` per exposure × outcome type. |
| `Meta_Sumstats/<analysis>/<analysis>.phewas_meta.plot_input.tsv.gz` | `variant_id` (the exposure), `phenotype`, `p_value` for the PheWAS plots. |
| `Plots/<analysis>.phewas_plot.<exposure>.png` | PheWAS (Manhattan-style) plot of the meta p-values per exposure; only when `pheno_descriptions_file` is set. `<analysis>.Plotting.log.txt` is the plotting log. |
| `Plots/<analysis>.forest.<exposure>.<outcome>.png` | Forest plot for each of the top `forest_plot_max` hits. |
| `Summary/phewas_meta_all_suggestive.csv` | All analyses' top hits in one table. |
| `Summary/phewas_meta_test_counts.csv`, `Summary/phewas_meta_sample_sizes.csv` | All analyses' test counts and sample sizes. |
| `Summary/phewas_meta_harmonization_log.csv` | All harmonization log rows. |
| `Summary/phewas_meta_cohorts.csv` | `meta_cohort, biobank, run_cohort`: one row per biobank × run cohort used by an analysis. |
| `Summary/saige_phewas_meta_params.json` | The run's parameters. |
| `PheWAS_Meta_Report.zip` | HTML report (see [Report](#report)); only with `pheno_descriptions_file`. |

## Report

`PheWAS_Meta_Report.zip` (in the launch directory) holds a self-contained HTML report: unzip it and open `PheWAS_Meta_Report/index.html`. It is made only when `pheno_descriptions_file` is set (as in Single-Variant PheWAS), because it includes the PheWAS plots. Sections, in order:

1. **Run Parameters**: the run's parameters (from `Summary/saige_phewas_meta_params.json`).
2. **Analyses**: the meta cohorts of each analysis, with the biobank and run cohort each came from. Biobank and run cohort come from `Summary/phewas_meta_cohorts.csv`.
3. **Harmonization Log**: every renamed, dropped or missing item. `missing_file` and `empty_file` rows are highlighted because those exposures are absent from the meta-analysis.
4. **Top Hits**: `phewas_meta_all_suggestive.csv` (paginated).
5. **PheWAS Plots** and **Forest Plots**: the PNGs in `Plots/`, grouped by analysis; click an image to enlarge it.

The zip also contains `Plots/*.png` and the four Summary CSVs.

## Plots

- **PheWAS plots** reuse the Single-Variant PheWAS plotting script on `<analysis>.phewas_meta.plot_input.tsv.gz` (`p_value` is the `meta_p_col` p-value). They are made only when `pheno_descriptions_file` is set; `varphewas_plot_max_labels` and `varphewas_plot_max_labels_per_category` control labelling.
- **Forest plots** are drawn for the first `forest_plot_max` rows of the filtered file (already sorted by `meta_p_col`). Each has one row per cohort with an estimate (beta with a 95% CI of 1.96 x SE; cohorts without a result are omitted) and an "IVW meta" diamond. The title gives `p_ivw`, `p_stouffer`, `p_het` and I². File names use the exposure and outcome with characters other than letters, numbers, `.`, `_` and `-` replaced by `_`; if two hits map to the same name (e.g. outcomes `A/B` and `A_B`), the later one gets an 8-character hash suffix (`<analysis>.forest.<exposure>.<outcome>.<hash>.png`). An analysis with an empty filtered file gets no forest plots. Forest plots always run.
- Both plot steps are `safe_to_skip`: with `skipPostprocessingErrors = true` a plotting failure does not stop the run.

# Caveats

- **Overlapping meta groups are not independent.** You define the analyses; if two analyses share cohorts (for example `ALL` and `AFR`), or a biobank's run cohorts overlap in samples (for example `ALL` and `AFR` inside one biobank listed as two meta cohorts of one analysis), the meta-analysis double-counts those samples and its p-values are too small. Put each sample in at most one meta cohort per analysis.
- Fixed-effect IVW assumes one true effect. Check `q_het` / `p_het` / `i2_het` before interpreting a pooled beta.
- The exposure-mode caveats of each SVP run (variance ratio calibration and approximate effect sizes; see [Caveats](Single-Variant_PheWAS_docs.md#caveats)) carry over to the pooled estimates.
- Only chromosome-1 exposure-mode sumstats (`*.chr1.variant_phewas.saige.gz`, where SVP writes exposure-mode results) are read; variant-mode SVP meta-analysis is not supported.

# Testing

`unit_tests/run_phewas_meta.sh` runs four cases against the synthetic two-biobank fixture in `test_data/PheWAS_META` (accepts `--profile`, `--keep-work`, `--outdir`, `--extra-config`; also part of `unit_tests/run_all.sh`):

| Test | Expects |
|------|---------|
| `phewas_meta_two_biobanks` | All outputs exist; the `ALL` meta statistics match `test_data/PheWAS_META/expected_meta.tsv` within 1e-6; the harmonization log shows renames, the missing file, the empty file, self-test and not-in-study drops; the single-cohort `AFR_only` analysis has no top hits. |
| `phewas_meta_transform_mismatch` | Fails: BB2's run used `inverse_normal` but the study says `none` (`cont_exposure_transform mismatch` in the log). |
| `phewas_meta_bad_label` | Fails before any task: analysis `ALL` names the undefined label `BB3_X`. |
| `phewas_meta_missing_cohort` | Fails before any task: BB2's YAML maps the nonexistent run cohort `EURR`, which matches no sumstats files (`allow_missing_cohorts = false`). |
