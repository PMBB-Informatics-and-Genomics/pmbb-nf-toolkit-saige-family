import json
import zipfile
import argparse
import shutil
import pandas as pd
from pathlib import Path


CSS = """
body {
    font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
    margin: 0;
    padding: 0;
    box-sizing: border-box;
    background-color: #f5f7fa;
    color: #2d3748;
    line-height: 1.6;
}

.container {
    max-width: 1400px;
    margin: 0 auto;
    padding: 2rem;
    margin-left: 260px;
}

h1 {
    font-size: 2.5rem;
    font-weight: 700;
    color: #1a202c;
    margin-bottom: 2rem;
    padding-bottom: 1rem;
    border-bottom: 3px solid #3182ce;
}

h2 {
    font-size: 2rem;
    font-weight: 600;
    color: #2d3748;
    margin: 2rem 0;
    text-align: center;
}

h3 {
    font-size: 1.5rem;
    color: #2d3748;
    margin: 1.5rem 0;
    padding-bottom: 0.5rem;
    border-bottom: 2px solid #e2e8f0;
}

h4 {
    font-size: 1.2rem;
    color: #4a5568;
    margin: 1rem 0 0.5rem 0;
}

.plot-container {
    background: white;
    padding: 2rem;
    border-radius: 12px;
    box-shadow: 0 4px 6px rgba(0, 0, 0, 0.05);
    margin: 2rem 0;
}

.plot-controls {
    display: flex;
    justify-content: center;
    gap: 1rem;
    margin: 1rem 0;
    flex-wrap: wrap;
}

.plot-toggle {
    padding: 0.75rem 1.5rem;
    border: none;
    border-radius: 8px;
    font-size: 1rem;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.2s ease;
    background-color: #e2e8f0;
    color: #2d3748;
}

.plot-toggle:hover {
    transform: translateY(-1px);
    box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
}

.plot-toggle.active {
    background-color: #3182ce;
    color: white;
}

.plot-wrapper {
    display: flex;
    justify-content: center;
    align-items: center;
    margin: 1.5rem 0;
}

.plot-image {
    max-width: 800px;
    width: 100%;
    height: auto;
    border-radius: 8px;
    box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    transition: all 0.3s ease;
    cursor: pointer;
}

.plot-image.enlarged {
    max-width: 1200px;
}

.table-container {
    background: white;
    padding: 1.5rem;
    border-radius: 12px;
    box-shadow: 0 4px 6px rgba(0, 0, 0, 0.05);
    margin: 2rem 0;
    overflow-x: auto;
}

table {
    width: 100%;
    border-collapse: separate;
    border-spacing: 0;
    margin: 1rem 0;
}

th {
    background-color: #2d3748;
    color: white;
    font-weight: 600;
    padding: 1rem;
    text-align: left;
    position: sticky;
    top: 0;
}

td {
    padding: 1rem;
    border-bottom: 1px solid #e2e8f0;
    font-size: 0.95rem;
}

tr:hover {
    background-color: #f7fafc;
}

.pagination-controls {
    margin-top: 1rem;
    display: flex;
    justify-content: center;
    align-items: center;
    gap: 1rem;
}

.pagination-button {
    padding: 0.75rem 1.5rem;
    border: none;
    border-radius: 8px;
    background-color: #3182ce;
    color: white;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.2s ease;
}

.pagination-button:hover:not(:disabled) {
    background-color: #2c5282;
    transform: translateY(-1px);
}

.pagination-button:disabled {
    opacity: 0.5;
    cursor: not-allowed;
}

#page-info {
    font-weight: 500;
    color: #2d3748;
}

.side-menu {
    position: fixed;
    left: 0;
    top: 0;
    width: 250px;
    height: 100vh;
    background-color: #2d3748;
    padding: 2rem 0;
    color: white;
    overflow-y: auto;
    box-shadow: 2px 0 5px rgba(0, 0, 0, 0.1);
    z-index: 1000;
}

.side-menu a {
    display: block;
    padding: 1rem 1.5rem;
    color: #e2e8f0;
    text-decoration: none;
    transition: all 0.2s ease;
    font-weight: 500;
}

.side-menu a:hover {
    background-color: #4a5568;
    padding-left: 2rem;
}

.side-menu .submenu {
    background-color: #1a202c;
    overflow: hidden;
}

.cohort-group {
    border-bottom: 1px solid #4a5568;
}

.cohort-group > a {
    background-image: url("data:image/svg+xml;charset=US-ASCII,%3Csvg%20xmlns%3D%22http://www.w3.org/2000/svg%22%20width%3D%22292.4%22%20height%3D%22292.4%22%3E%3Cpath%20fill%3D%22%23ffffff%22%20d%3D%22M287%2069.4a17.6%2017.6%200%200%200-13-5.4H18.4c-5%200-9.3%201.8-12.9%205.4A17.6%2017.6%200%200%200%200%2082.2c0%205%201.8%209.3%205.4%2012.9l128%20127.9c3.6%203.6%207.8%205.4%2012.8%205.4s9.2-1.8%2012.8-5.4L287%2095c3.5-3.5%205.4-7.8%205.4-12.8%200-5-1.9-9.2-5.4-12.8z%22/%3E%3C/svg%3E");
    background-repeat: no-repeat;
    background-position: right 1rem center;
    background-size: 0.8em auto;
}

.subcohorts {
    background-color: #1a202c;
}

.subcohort-link {
    padding-left: 2.5rem !important;
    font-size: 0.9rem;
    color: #cbd5e0 !important;
}

.subcohort-link:hover {
    background-color: #4a5568;
    padding-left: 3rem !important;
}

.menu-toggle {
    display: none;
    position: fixed;
    top: 1rem;
    left: 1rem;
    z-index: 1001;
    padding: 0.5rem;
    background-color: #2d3748;
    border-radius: 4px;
    cursor: pointer;
}

.menu-toggle span {
    display: block;
    width: 25px;
    height: 3px;
    background-color: white;
    margin: 5px 0;
    transition: all 0.3s ease;
}

@media (max-width: 1024px) {
    .side-menu {
        transform: translateX(-100%);
        transition: transform 0.3s ease;
    }

    .side-menu.active {
        transform: translateX(0);
    }

    .menu-toggle {
        display: block;
    }

    .container {
        margin-left: 0;
        padding: 1rem;
    }

    .plot-image {
        max-width: 100%;
    }

    .plot-image.enlarged {
        max-width: 100%;
    }
}
"""

JS = """
function toggleMenu() {
    document.querySelector('.side-menu').classList.toggle('active');
}

function toggleSubmenu(id) {
    const el = document.getElementById(id);
    el.style.display = (el.style.display === 'none') ? 'block' : 'none';
}

function toggleSubCohorts(pheno) {
    const el = document.getElementById(`${pheno}-subcohorts`);
    el.style.display = (el.style.display === 'none') ? 'block' : 'none';
}

// ---- Plot type toggle (manhattan / qq) ----
function togglePlotType(sectionId, plotType) {
    const section = document.getElementById(sectionId);
    if (!section) return;
    section.querySelectorAll('.plot-toggle').forEach(btn => {
        btn.classList.toggle('active', btn.dataset.plotType === plotType);
    });
    section.querySelectorAll('.plot-wrapper').forEach(wrapper => {
        wrapper.style.display = (wrapper.dataset.plotType === plotType) ? 'flex' : 'none';
    });
}

// ---- Table pagination ----
let currentPage = 1;
let rowsPerPage = 10;

function showCurrentPage() {
    const startIdx = (currentPage - 1) * rowsPerPage;
    const endIdx = startIdx + rowsPerPage;
    const table = document.getElementById('results-table');
    if (!table) return;
    const allRows = Array.from(table.querySelector('tbody').getElementsByTagName('tr'));
    allRows.forEach(row => row.style.display = 'none');
    allRows.slice(startIdx, endIdx).forEach(row => { row.style.display = ''; });
    updatePaginationInfo(allRows.length);
}

function updatePaginationInfo(total) {
    const totalPages = Math.ceil(total / rowsPerPage);
    const pageInfo = document.getElementById('page-info');
    const prev = document.getElementById('prev-page');
    const next = document.getElementById('next-page');
    if (pageInfo) pageInfo.textContent = `Page ${currentPage} of ${totalPages}`;
    if (prev) prev.disabled = (currentPage <= 1);
    if (next) next.disabled = (currentPage >= totalPages);
}

function previousPage() {
    if (currentPage > 1) { currentPage--; showCurrentPage(); }
}

function nextPage() {
    const table = document.getElementById('results-table');
    if (!table) return;
    const totalRows = table.querySelector('tbody').getElementsByTagName('tr').length;
    if (currentPage < Math.ceil(totalRows / rowsPerPage)) { currentPage++; showCurrentPage(); }
}

function toggleImageSize(img) {
    img.classList.toggle('enlarged');
}

document.addEventListener('DOMContentLoaded', () => { showCurrentPage(); });
"""


class SaigeGwasReportGenerator:
    def __init__(self, manifest_path, output_zip):
        with open(manifest_path) as f:
            self.manifest = json.load(f)

        self.output_zip = Path(output_zip)
        self.report_name = self.output_zip.stem
        self.output_dir = Path(self.report_name)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / 'Plots').mkdir(exist_ok=True)

        self.cohort_list = self.manifest['cohort_list']
        self.bin_pheno_list = self.manifest.get('bin_pheno_list', [])
        self.quant_pheno_list = self.manifest.get('quant_pheno_list', [])
        self.survival_pheno_list = self.manifest.get('survival_pheno_list', [])
        self.all_phenos = self.bin_pheno_list + self.quant_pheno_list + self.survival_pheno_list

        self.gwas_df = pd.read_csv(self.manifest['gwas_summary_csv'])
        self.pheno_summaries_df = pd.read_csv(self.manifest['pheno_summaries_csv'])

    def _plot_rel_path(self, src_path):
        return f"Plots/{Path(src_path).name}"

    def _write_css(self):
        (self.output_dir / 'styles.css').write_text(CSS)

    def _create_sidebar(self):
        s = '<div class="side-menu">\n'
        s += '  <a href="index.html">Home</a>\n'
        s += '  <a href="phenotype_summary.html">Phenotype Summary</a>\n'
        s += "  <a href=\"#\" onclick=\"toggleSubmenu('results-filter-submenu')\">Results Filter</a>\n"
        s += '  <div id="results-filter-submenu" class="submenu">\n'
        for pheno in self.all_phenos:
            s += '    <div class="cohort-group">\n'
            s += f"      <a href=\"#\" onclick=\"toggleSubCohorts('{pheno}')\">{pheno}</a>\n"
            s += f'      <div id="{pheno}-subcohorts" class="subcohorts" style="display: none;">\n'
            for cohort in self.cohort_list:
                s += f'        <a href="{pheno}.{cohort}.gwas.html" class="subcohort-link">{cohort}</a>\n'
            s += '      </div>\n'
            s += '    </div>\n'
        s += '  </div>\n'
        s += '  <a href="method_summary.html">Analysis Logs</a>\n'
        s += '</div>\n'
        return s

    def _page_template(self, content, title="SAIGE GWAS Results Report"):
        return (
            '<!DOCTYPE html>\n<html lang="en">\n<head>\n'
            '    <meta charset="UTF-8">\n'
            '    <meta name="viewport" content="width=device-width, initial-scale=1.0">\n'
            f'    <title>{title}</title>\n'
            '    <link rel="stylesheet" href="styles.css">\n'
            '</head>\n<body>\n'
            '    <div class="menu-toggle" onclick="toggleMenu()">'
            '<span></span><span></span><span></span></div>\n'
            + self._create_sidebar()
            + '    <div class="container">\n'
            f'        <h1>{title}</h1>\n'
            + content
            + '    </div>\n'
            '<script>\n' + JS + '\n</script>\n'
            '</body>\n</html>\n'
        )

    def _df_to_html_table(self, df, table_id='results-table'):
        headers = df.columns.tolist()
        header_html = ''.join(f'<th>{h}</th>' for h in headers)
        rows_html = ''.join(
            '<tr>' + ''.join(f'<td>{row[c]}</td>' for c in headers) + '</tr>\n'
            for _, row in df.iterrows()
        )
        return (
            f'<div class="table-container"><table id="{table_id}">'
            f'<thead><tr>{header_html}</tr></thead>'
            f'<tbody>{rows_html}</tbody>'
            '</table></div>\n'
        )

    def _pagination_controls(self):
        return (
            '<div class="pagination-controls">'
            '<button id="prev-page" class="pagination-button" onclick="previousPage()">Previous</button>'
            '<span id="page-info">Page 1</span>'
            '<button id="next-page" class="pagination-button" onclick="nextPage()">Next</button>'
            '</div>\n'
        )

    def _filter_hits(self, df, cohort, pheno):
        if 'COHORT' in df.columns and 'PHENO' in df.columns:
            return df[(df['COHORT'] == cohort) & (df['PHENO'] == pheno)]
        return pd.DataFrame()

    def generate_index_page(self):
        content = """
        <div id="default-view">
            <p style="text-align: left;">
                This report presents findings from a SAIGE GWAS (Genome-Wide Association Study) analysis.
                Single-variant association tests were performed across the specified phenotypes and cohorts.
            </p>
            <h2 style="text-align: left;">Getting Started</h2>
            <ul style="text-align: left;">
                <li>Use the <b>Results Filter</b> in the sidebar to navigate to a phenotype and cohort result.</li>
                <li>Check the <b>Phenotype Summary</b> for trait distributions and sample counts.</li>
                <li>Each result page shows Manhattan and QQ plots alongside top hits above the suggestive threshold.</li>
                <li>Refer to <b>Analysis Logs</b> for the full pipeline parameter set.</li>
            </ul>
        </div>
        """
        (self.output_dir / 'index.html').write_text(
            self._page_template(content, 'SAIGE GWAS Results Report')
        )
        print("  index.html")

    def generate_phenotype_summary(self):
        table_html = self._df_to_html_table(self.pheno_summaries_df, table_id='pheno-summary-table')

        plot_sections = ''
        for pheno, plots in self.manifest.get('pheno_summary_plots', {}).items():
            imgs = ''.join(
                f'<div class="plot-wrapper">'
                f'<img src="{self._plot_rel_path(p)}" alt="{pheno} summary" '
                f'class="plot-image" onclick="toggleImageSize(this)"></div>\n'
                for p in plots
            )
            plot_sections += (
                f'<div class="plot-container"><h3>{pheno}</h3>{imgs}</div>\n'
            )

        content = (
            '<h2>Phenotype Summary Statistics</h2>\n'
            + table_html
            + '<h2>Phenotype Distribution Plots</h2>\n'
            + (plot_sections if plot_sections else '<p>No phenotype summary plots found.</p>')
        )
        (self.output_dir / 'phenotype_summary.html').write_text(
            self._page_template(content, 'Phenotype Summary')
        )
        print("  phenotype_summary.html")

    def generate_gwas_page(self, cohort, pheno):
        plots = self.manifest['gwas_plots'].get(cohort, {}).get(pheno, {})
        section_id = f'gwas-{pheno}-{cohort}'

        plot_content = ''
        if plots:
            manhattan_src = plots.get('manhattan')
            qq_src = plots.get('qq')
            plot_content = (
                f'<div class="plot-container" id="{section_id}">\n'
                f'<h3>GWAS Plots — {pheno} in {cohort}</h3>\n'
                '<div class="plot-controls">'
                f'<button class="plot-toggle active" data-plot-type="manhattan" '
                f'onclick="togglePlotType(\'{section_id}\', \'manhattan\')">Manhattan Plot</button>'
                f'<button class="plot-toggle" data-plot-type="qq" '
                f'onclick="togglePlotType(\'{section_id}\', \'qq\')">QQ Plot</button>'
                '</div>\n'
            )
            for pt_type, src, display in [('manhattan', manhattan_src, 'flex'), ('qq', qq_src, 'none')]:
                if src:
                    plot_content += (
                        f'<div class="plot-wrapper" data-plot-type="{pt_type}" style="display: {display}">'
                        f'<img src="{self._plot_rel_path(src)}" alt="{cohort} {pheno} {pt_type}" '
                        f'class="plot-image" onclick="toggleImageSize(this)"></div>\n'
                    )
            plot_content += '</div>\n'

        df = self._filter_hits(self.gwas_df, cohort, pheno)
        if df.empty:
            table_html = '<div class="table-container"><p>No significant hits found above the p-value threshold.</p></div>\n'
            pagination = ''
        else:
            table_html = self._df_to_html_table(df)
            pagination = self._pagination_controls()

        content = (
            f'<h2>GWAS Results — {pheno} in {cohort}</h2>\n'
            + plot_content
            + '<h3>Top GWAS Hits</h3>\n'
            + table_html + pagination
        )
        (self.output_dir / f'{pheno}.{cohort}.gwas.html').write_text(
            self._page_template(content, f'GWAS — {pheno} in {cohort}')
        )
        print(f"  {pheno}.{cohort}.gwas.html")

    def generate_method_summary(self):
        params = self.manifest.get('params', {})

        def fmt(v):
            if isinstance(v, list):
                return ', '.join(str(x) for x in v) or '—'
            if isinstance(v, dict):
                return '<br>'.join(f'{k} → {val}' for k, val in v.items())
            return str(v) if v is not None else '—'

        def fmt_key(k):
            return k.replace('_', ' ').title()

        body = ''
        for section_key, section_data in params.items():
            if not isinstance(section_data, dict):
                continue
            rows = ''.join(
                f'<tr><td><b>{fmt_key(k)}</b></td><td>{fmt(v)}</td></tr>\n'
                for k, v in section_data.items()
            )
            body += (
                f'<h3>{fmt_key(section_key)}</h3>'
                '<div class="table-container"><table>'
                '<thead><tr><th>Parameter</th><th>Value</th></tr></thead>'
                f'<tbody>{rows}</tbody>'
                '</table></div>\n'
            )

        content = (
            '<div id="method-summary">'
            '<h2>SAIGE GWAS Methods Summary</h2>'
            + body
            + '</div>\n'
        )
        (self.output_dir / 'method_summary.html').write_text(
            self._page_template(content, 'Analysis Logs')
        )
        print("  method_summary.html")

    def _copy_assets(self):
        all_plot_sources = set()
        for plots in self.manifest.get('pheno_summary_plots', {}).values():
            all_plot_sources.update(plots)
        for cohort_plots in self.manifest.get('gwas_plots', {}).values():
            for pheno_plots in cohort_plots.values():
                all_plot_sources.update(pheno_plots.values())

        for src in all_plot_sources:
            src_path = Path(src)
            if src_path.exists():
                shutil.copy(src_path, self.output_dir / 'Plots' / src_path.name)
            else:
                print(f"  Warning: plot not found: {src}")

        shutil.copy(self.manifest['gwas_summary_csv'],
                    self.output_dir / 'saige_gwas_suggestive.csv')
        shutil.copy(self.manifest['pheno_summaries_csv'],
                    self.output_dir / 'pheno_summaries.csv')

    def _create_zip(self):
        with zipfile.ZipFile(self.output_zip, 'w', zipfile.ZIP_DEFLATED) as zf:
            for file_path in self.output_dir.rglob('*'):
                if file_path.is_file():
                    zf.write(file_path, file_path.relative_to(self.output_dir.parent))
        print(f"Report written to {self.output_zip}")

    def generate_all(self):
        print("Generating SAIGE GWAS report...")
        self._write_css()
        self.generate_index_page()
        self.generate_phenotype_summary()
        for pheno in self.all_phenos:
            for cohort in self.cohort_list:
                self.generate_gwas_page(cohort, pheno)
        self.generate_method_summary()
        self._copy_assets()
        self._create_zip()
        print("Done.")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a portable SAIGE GWAS HTML report zip from a manifest file."
    )
    parser.add_argument('--manifest', required=True,
                        help='Path to saige_gwas_manifest.json')
    parser.add_argument('--output_zip', default='SAIGE_GWAS_Report.zip',
                        help='Output zip file path (default: SAIGE_GWAS_Report.zip)')
    args = parser.parse_args()

    generator = SaigeGwasReportGenerator(args.manifest, args.output_zip)
    generator.generate_all()


if __name__ == '__main__':
    main()
