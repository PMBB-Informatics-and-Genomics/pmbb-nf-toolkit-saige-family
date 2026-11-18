import json
import html
import shutil
import zipfile
import argparse
from pathlib import Path

import pandas as pd


CSS = """
body { font-family: 'Segoe UI', system-ui, -apple-system, sans-serif; margin: 0; background-color: #f5f7fa;
       color: #2d3748; line-height: 1.6; }
.container { max-width: 1400px; margin: 0 auto; padding: 2rem; margin-left: 260px; }
h1 { font-size: 2.2rem; font-weight: 700; color: #1a202c; margin-bottom: 2rem; padding-bottom: 1rem;
     border-bottom: 3px solid #3182ce; }
h2 { font-size: 1.7rem; font-weight: 600; color: #2d3748; margin: 3rem 0 1rem; }
h3 { font-size: 1.3rem; color: #2d3748; margin: 1.5rem 0; padding-bottom: 0.5rem; border-bottom: 2px solid #e2e8f0; }
.plot-container { background: white; padding: 1.5rem; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.05);
                  margin: 1.5rem 0; }
.plot-grid { display: flex; flex-wrap: wrap; gap: 1.5rem; justify-content: center; }
.plot-wrapper { flex: 1 1 520px; max-width: 800px; text-align: center; }
.plot-image { max-width: 100%; height: auto; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); cursor: pointer; }
.plot-image.enlarged { max-width: 1200px; width: 100%; }
.table-container { background: white; padding: 1.5rem; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.05);
                   margin: 1rem 0; overflow-x: auto; }
table { width: 100%; border-collapse: separate; border-spacing: 0; }
th { background-color: #2d3748; color: white; font-weight: 600; padding: 0.6rem 0.8rem; text-align: left; }
td { padding: 0.5rem 0.8rem; border-bottom: 1px solid #e2e8f0; font-size: 0.92rem; }
tr:hover { background-color: #f7fafc; }
tr.warn-row td { background-color: #fefcbf; }
.pagination-controls { margin-top: 0.8rem; display: flex; justify-content: center; align-items: center; gap: 1rem; }
.pagination-button { padding: 0.5rem 1.2rem; border: none; border-radius: 8px; background-color: #3182ce; color: white;
                     cursor: pointer; }
.pagination-button:disabled { opacity: 0.5; cursor: not-allowed; }
.side-menu { position: fixed; left: 0; top: 0; width: 250px; height: 100vh; background-color: #2d3748;
             padding: 2rem 0; overflow-y: auto; z-index: 1000; }
.side-menu a { display: block; padding: 0.8rem 1.5rem; color: #e2e8f0; text-decoration: none; font-weight: 500; }
.side-menu a:hover { background-color: #4a5568; }
@media (max-width: 1024px) {
    .side-menu { display: none; }
    .container { margin-left: 0; padding: 1rem; }
}
"""

JS = """
const ROWS_PER_PAGE = 15;
function paginate(table) {
    const rows = Array.from(table.querySelector('tbody').rows);
    const pages = Math.max(1, Math.ceil(rows.length / ROWS_PER_PAGE));
    if (rows.length <= ROWS_PER_PAGE) return;
    let page = 1;
    const bar = document.createElement('div');
    bar.className = 'pagination-controls';
    bar.innerHTML = '<button class="pagination-button">Previous</button><span></span>' +
                    '<button class="pagination-button">Next</button>';
    table.closest('.table-container').after(bar);
    const [prev, info, next] = [bar.children[0], bar.children[1], bar.children[2]];
    function show() {
        rows.forEach((r, i) => { r.style.display = (i >= (page - 1) * ROWS_PER_PAGE && i < page * ROWS_PER_PAGE) ? '' : 'none'; });
        info.textContent = 'Page ' + page + ' of ' + pages + ' (' + rows.length + ' rows)';
        prev.disabled = page <= 1; next.disabled = page >= pages;
    }
    prev.onclick = () => { page--; show(); };
    next.onclick = () => { page++; show(); };
    show();
}
function toggleImageSize(img) { img.classList.toggle('enlarged'); }
document.addEventListener('DOMContentLoaded', () => { document.querySelectorAll('table.paginated').forEach(paginate); });
"""

# harmonization log actions that mean data is absent: highlighted in the log table
WARN_ACTIONS = ('missing_file', 'empty_file')


def esc(value):
    return html.escape(str(value), quote=True)


def read_csv_or_empty(path):
    try:
        return pd.read_csv(path, dtype=str, keep_default_na=False)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


def fmt_cell(value):
    """Floats get 6 significant digits; everything else is shown as is."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if value == '' or str(value).lstrip('+-').isdigit():
        return str(value)
    return f'{f:.6g}'


# identifier columns of the top hits table: never reformatted as numbers (outcome names such as 250.2)
ID_COLS = ('analysis', 'exposure', 'outcome', 'exposure_type', 'outcome_type', 'direction')


def df_to_html_table(df, table_id, paginated=False, warn_col=None, numeric=False):
    header = ''.join(f'<th>{esc(c)}</th>' for c in df.columns)
    rows = []
    for _, row in df.iterrows():
        warn = warn_col is not None and row[warn_col] in WARN_ACTIONS
        cells = ''.join(f'<td>{esc(fmt_cell(row[c]) if numeric and c not in ID_COLS else row[c])}</td>' for c in df.columns)
        rows.append(f'<tr class="warn-row">{cells}</tr>\n' if warn else f'<tr>{cells}</tr>\n')
    cls = ' class="paginated"' if paginated else ''
    return (f'<div class="table-container"><table id="{esc(table_id)}"{cls}><thead><tr>{header}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>\n')


def fmt_param(value):
    if value is None or value == '':
        return '-'
    if isinstance(value, list):
        return ', '.join(str(v) for v in value) or '-'
    if isinstance(value, dict):
        return '; '.join(f'{k}: {fmt_param(v)}' for k, v in value.items()) or '-'
    return str(value)


def render_params(params):
    rows = ''.join(f'<tr><td><b>{esc(k)}</b></td><td>{esc(fmt_param(v))}</td></tr>\n' for k, v in params.items())
    return ('<div class="table-container"><table><thead><tr><th>Parameter</th><th>Value</th></tr></thead>'
            f'<tbody>{rows}</tbody></table></div>\n')


def render_analyses(analyses):
    rows = []
    for name, cohorts in analyses.items():
        for c in cohorts:
            rows.append({'analysis': name, 'meta cohort': c['cohort'], 'biobank': c.get('biobank') or '?',
                         'run cohort': c.get('run_cohort') or '?'})
    if not rows:
        return '<p>No analyses.</p>\n'
    return df_to_html_table(pd.DataFrame(rows), 'analyses-table')


def render_plot_gallery(plots_by_analysis, label):
    if not plots_by_analysis:
        return f'<p>No {esc(label)} were produced.</p>\n'
    out = ''
    for analysis, pngs in plots_by_analysis.items():
        imgs = ''.join(
            f'<div class="plot-wrapper"><img src="Plots/{esc(Path(p).name)}" alt="{esc(Path(p).stem)}" '
            f'class="plot-image" onclick="toggleImageSize(this)"></div>\n'
            for p in pngs)
        out += f'<div class="plot-container"><h3>{esc(analysis)}</h3>\n<div class="plot-grid">{imgs}</div></div>\n'
    return out


def build_report_html(manifest, title='PheWAS-META Report'):
    hits = read_csv_or_empty(manifest['top_hits'])
    log = read_csv_or_empty(manifest['harmonization_log'])

    hits_html = (df_to_html_table(hits, 'top-hits-table', paginated=True, numeric=True)
                 if len(hits) else '<p>No top hits above the p-value cutoff.</p>\n')
    if len(log) and 'action' in log.columns:
        log_html = df_to_html_table(log, 'harmonization-log-table', paginated=True, warn_col='action')
        n_warn = int(log['action'].isin(WARN_ACTIONS).sum())
        log_note = (f'<p>{n_warn} row(s) with <code>missing_file</code> or <code>empty_file</code> are highlighted: '
                    'those exposure files were absent or empty and are not part of the meta-analysis.</p>\n')
    else:
        log_html, log_note = '<p>The harmonization log is empty: nothing was renamed, dropped or missing.</p>\n', ''

    sections = [
        ('params', 'Run Parameters', render_params(manifest.get('params', {}))),
        ('analyses', 'Analyses',
         '<p>Meta cohorts in each analysis, with the biobank and run cohort each came from.</p>\n' + render_analyses(manifest.get('analyses', {}))),
        ('harmonization', 'Harmonization Log', log_note + log_html),
        ('hits', 'Top Hits', hits_html),
        ('phewas', 'PheWAS Plots', render_plot_gallery(manifest.get('phewas_plots', {}), 'PheWAS plots')),
        ('forest', 'Forest Plots', render_plot_gallery(manifest.get('forest_plots', {}), 'forest plots')),
    ]
    menu = ''.join(f'  <a href="#{sid}">{esc(h)}</a>\n' for sid, h, _ in sections)
    body = ''.join(f'<h2 id="{sid}">{esc(h)}</h2>\n{content}' for sid, h, content in sections)
    return (
        '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="UTF-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1.0">\n'
        f'<title>{esc(title)}</title>\n<style>{CSS}</style>\n</head>\n<body>\n'
        f'<div class="side-menu">\n{menu}</div>\n'
        f'<div class="container">\n<h1>{esc(title)}</h1>\n{body}</div>\n'
        f'<script>{JS}</script>\n</body>\n</html>\n'
    )


def write_html_and_zip(manifest, output_zip):
    output_zip = Path(output_zip)
    out_dir = Path(output_zip.stem)
    (out_dir / 'Plots').mkdir(parents=True, exist_ok=True)

    (out_dir / 'index.html').write_text(build_report_html(manifest))
    for group in ('phewas_plots', 'forest_plots'):
        for pngs in manifest.get(group, {}).values():
            for p in pngs:
                if Path(p).exists():
                    shutil.copy(p, out_dir / 'Plots' / Path(p).name)
                else:
                    print(f"  Warning: plot not found: {p}")
    for key in ('top_hits', 'test_counts', 'sample_sizes', 'harmonization_log'):
        if manifest.get(key) and Path(manifest[key]).exists():
            shutil.copy(manifest[key], out_dir / Path(manifest[key]).name)

    with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(out_dir.rglob('*')):
            if f.is_file():
                zf.write(f, f.relative_to(out_dir.parent))
    print(f"Report written to {output_zip}")


def main():
    parser = argparse.ArgumentParser(description="Generate a portable PheWAS-META HTML report zip.")
    parser.add_argument('--manifest', required=True, help='Path to phewas_meta_manifest.json')
    parser.add_argument('--output_zip', default='PheWAS_Meta_Report.zip', help='Output zip file path')
    args = parser.parse_args()
    with open(args.manifest) as f:
        manifest = json.load(f)
    write_html_and_zip(manifest, args.output_zip)


if __name__ == '__main__':
    main()
