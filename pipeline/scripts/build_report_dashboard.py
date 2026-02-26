import csv
import json
import pathlib
import sys


def load_tsv(path):
    with pathlib.Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def parse_num(value):
    try:
        return float(value)
    except Exception:
        return None


def render_dashboard_html(metrics, signals):
    data_blob = json.dumps({'metrics': metrics, 'signals': signals})
    return f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>EBV Report Dashboard</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 20px; color: #1f2937; }}
    h1 {{ margin: 0 0 10px 0; font-size: 24px; }}
    .grid {{ display: grid; grid-template-columns: repeat(4, minmax(180px, 1fr)); gap: 10px; margin: 14px 0 20px 0; }}
    .card {{ border: 1px solid #d1d5db; border-radius: 8px; padding: 10px; background: #f9fafb; }}
    .label {{ font-size: 12px; color: #4b5563; margin-bottom: 6px; }}
    .value {{ font-size: 18px; font-weight: 600; }}
    .controls {{ display: flex; gap: 10px; margin-bottom: 10px; }}
    select, input {{ padding: 6px 8px; border: 1px solid #d1d5db; border-radius: 6px; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
    th, td {{ border: 1px solid #e5e7eb; padding: 6px 8px; text-align: left; }}
    th {{ background: #f3f4f6; cursor: pointer; }}
    .muted {{ color: #6b7280; font-size: 12px; }}
  </style>
</head>
<body>
  <h1>EBV Multi-Omics Dashboard</h1>
  <p class=\"muted\">Phase 6 Wave 2 interactive package generated from report bundle outputs.</p>

  <div class=\"grid\" id=\"metric-cards\"></div>

  <div class=\"controls\">
    <select id=\"sourceFilter\">
      <option value=\"all\">All sources</option>
      <option value=\"wgs_assoc\">WGS association</option>
      <option value=\"rna_assoc\">RNA association</option>
    </select>
    <input id=\"featureSearch\" placeholder=\"Search feature/gene\" />
  </div>

  <table>
    <thead>
      <tr>
        <th data-key=\"source\">source</th>
        <th data-key=\"pathway\">pathway</th>
        <th data-key=\"layer\">layer</th>
        <th data-key=\"feature_id\">feature_id</th>
        <th data-key=\"gene_id\">gene_id</th>
        <th data-key=\"effect_size\">effect_size</th>
        <th data-key=\"p_value\">p_value</th>
        <th data-key=\"p_value_fdr\">p_value_fdr</th>
        <th data-key=\"n_samples\">n_samples</th>
      </tr>
    </thead>
    <tbody id=\"signalsBody\"></tbody>
  </table>

<script>
const DATA = {data_blob};
let sortKey = 'p_value';
let sortAsc = true;

function metricCards() {{
  const cards = [
    ['WGS EBV coverage', DATA.metrics?.wgs_assoc?.ebv_model_coverage_rate],
    ['WGS host coverage', DATA.metrics?.wgs_assoc?.host_model_coverage_rate],
    ['RNA EBV genes', DATA.metrics?.rna_assoc?.ebv_gene_count],
    ['Merged top signals', DATA.metrics?.bundle?.merged_top_signal_count],
  ];
  const root = document.getElementById('metric-cards');
  root.innerHTML = '';
  for (const [label, value] of cards) {{
    const card = document.createElement('div');
    card.className = 'card';
    card.innerHTML = `<div class=\"label\">${{label}}</div><div class=\"value\">${{value ?? 'NA'}}</div>`;
    root.appendChild(card);
  }}
}}

function numeric(v) {{
  const x = parseFloat(v);
  return Number.isFinite(x) ? x : null;
}}

function filteredSignals() {{
  const source = document.getElementById('sourceFilter').value;
  const search = document.getElementById('featureSearch').value.toLowerCase().trim();
  return (DATA.signals || []).filter(row => {{
    if (source !== 'all' && row.source !== source) return false;
    if (!search) return true;
    return `${{row.feature_id || ''}} ${{row.gene_id || ''}}`.toLowerCase().includes(search);
  }});
}}

function sortedSignals(rows) {{
  return rows.slice().sort((a, b) => {{
    const ax = numeric(a[sortKey]);
    const bx = numeric(b[sortKey]);
    let cmp = 0;
    if (ax !== null && bx !== null) cmp = ax - bx;
    else cmp = String(a[sortKey] || '').localeCompare(String(b[sortKey] || ''));
    return sortAsc ? cmp : -cmp;
  }});
}}

function renderTable() {{
  const body = document.getElementById('signalsBody');
  body.innerHTML = '';
  const rows = sortedSignals(filteredSignals());
  for (const row of rows) {{
    const tr = document.createElement('tr');
    tr.innerHTML = `
      <td>${{row.source || ''}}</td>
      <td>${{row.pathway || ''}}</td>
      <td>${{row.layer || ''}}</td>
      <td>${{row.feature_id || ''}}</td>
      <td>${{row.gene_id || ''}}</td>
      <td>${{row.effect_size || ''}}</td>
      <td>${{row.p_value || ''}}</td>
      <td>${{row.p_value_fdr || ''}}</td>
      <td>${{row.n_samples || ''}}</td>`;
    body.appendChild(tr);
  }}
}}

for (const th of document.querySelectorAll('th[data-key]')) {{
  th.addEventListener('click', () => {{
    const key = th.getAttribute('data-key');
    if (sortKey === key) sortAsc = !sortAsc;
    else {{ sortKey = key; sortAsc = true; }}
    renderTable();
  }});
}}

document.getElementById('sourceFilter').addEventListener('change', renderTable);
document.getElementById('featureSearch').addEventListener('input', renderTable);

metricCards();
renderTable();
</script>
</body>
</html>
"""


def main():
    metrics_path = pathlib.Path(sys.argv[1])
    top_signals_path = pathlib.Path(sys.argv[2])

    metrics = json.loads(metrics_path.read_text())
    signals = load_tsv(top_signals_path)

    dashboard_data = {
        'metrics': metrics,
        'signal_count': len(signals),
    }
    pathlib.Path('report_dashboard_data.json').write_text(json.dumps(dashboard_data, indent=2))
    pathlib.Path('report_dashboard.html').write_text(render_dashboard_html(metrics, signals))


if __name__ == '__main__':
    main()
