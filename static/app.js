const statsEl = document.getElementById('stats');
const resultsBody = document.getElementById('resultsBody');
const featureBody = document.getElementById('featureBody');
const detailMeta = document.getElementById('detailMeta');

const gseInput = document.getElementById('gseInput');
const gsmInput = document.getElementById('gsmInput');
const qInput = document.getElementById('qInput');
const limitInput = document.getElementById('limitInput');

let currentRows = [];
let pickedGsm = '';

function card(k, v) {
  return `<div class="stat-card"><div class="k">${k}</div><div class="v">${v}</div></div>`;
}

function meta(k, v) {
  return `<div class="meta-item"><div class="k">${k}</div><div class="v">${v || '-'}</div></div>`;
}

async function loadStats() {
  const resp = await fetch('/api/stats');
  const data = await resp.json();
  statsEl.innerHTML = [
    card('样本数', data.samples),
    card('GSE 数', data.gses),
    card('Gene/Value 行数', data.gene_like_rows),
    card('更新时间', data.db_updated_at),
  ].join('');
}

function renderRows(rows) {
  currentRows = rows;
  resultsBody.innerHTML = rows.map((r) => `
    <tr data-gsm="${r.gsm_id}" class="${pickedGsm === r.gsm_id ? 'pick' : ''}">
      <td>${r.gse_id || '-'}</td>
      <td>${r.gsm_id}</td>
      <td>${r.sample_title || '-'}</td>
      <td>${r.treatment || '-'}</td>
      <td>${r.genotype || '-'}</td>
      <td>${r.feature_count}</td>
      <td>${r.parse_quality || '-'}</td>
    </tr>
  `).join('');

  resultsBody.querySelectorAll('tr').forEach((tr) => {
    tr.addEventListener('click', () => {
      pickedGsm = tr.dataset.gsm;
      renderRows(currentRows);
      loadDetail(pickedGsm);
    });
  });
}

async function search() {
  const params = new URLSearchParams({
    gse: gseInput.value.trim(),
    gsm: gsmInput.value.trim(),
    q: qInput.value.trim(),
    limit: String(limitInput.value || 100),
  });

  const resp = await fetch(`/api/search?${params.toString()}`);
  const data = await resp.json();
  renderRows(data.rows || []);

  if ((data.rows || []).length > 0) {
    pickedGsm = data.rows[0].gsm_id;
    renderRows(data.rows);
    loadDetail(pickedGsm);
  } else {
    detailMeta.innerHTML = '';
    featureBody.innerHTML = '<tr><td colspan="3">没有匹配结果</td></tr>';
  }
}

async function loadDetail(gsm) {
  if (!gsm) return;
  const resp = await fetch(`/api/sample/${gsm}`);
  if (!resp.ok) {
    detailMeta.innerHTML = '<div>未找到该样本</div>';
    featureBody.innerHTML = '';
    return;
  }

  const data = await resp.json();
  const s = data.sample;
  detailMeta.innerHTML = [
    meta('GSE', s.gse_id),
    meta('GSM', s.gsm_id),
    meta('样本名', s.sample_title),
    meta('treatment', s.treatment),
    meta('genotype', s.genotype),
    meta('解析质量', s.parse_quality),
    meta('原始条件', s.raw_characteristics),
  ].join('');

  const rows = data.gene_like_rows || [];
  featureBody.innerHTML = rows.length
    ? rows.map((x) => `<tr><td>${x.gene}</td><td>${x.value}</td><td>${x.source}</td></tr>`).join('')
    : '<tr><td colspan="3">暂无明细</td></tr>';
}

async function bootstrap() {
  await loadStats();
  await search();
}

document.getElementById('searchBtn').addEventListener('click', search);
document.getElementById('resetBtn').addEventListener('click', async () => {
  gseInput.value = '';
  gsmInput.value = '';
  qInput.value = '';
  limitInput.value = 100;
  pickedGsm = '';
  await search();
});

bootstrap();
