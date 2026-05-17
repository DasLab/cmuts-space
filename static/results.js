// Results page: per-group + per-sequence plot navigation.

(function () {
  const meta = JSON.parse(document.getElementById('meta-data').textContent);
  const jobId = window.location.pathname.split('/').filter(Boolean)[1];

  const groupSelect = document.getElementById('group-select');
  const seqSelect = document.getElementById('seq-select');
  const seqWrapper = document.getElementById('seq-select-wrapper');
  const statsTable = document.getElementById('stats-table');
  const tiles = Array.from(document.querySelectorAll('.plot-tile'));

  // Plots whose data is per-reference (sequence selector applies).
  const PER_REF = new Set(['profile', 'mi', 'correlation', 'pairwise_coverage']);

  function emptyTile(tile, msg) {
    const plot = tile.querySelector('.plot');
    plot.innerHTML = `<p class="muted">${msg}</p>`;
  }

  function normalizeLayout(layout) {
    // Strip fixed width/height so the container drives sizing. Everything
    // else (theme, fonts, axes, colors) is left exactly as cmuts emits.
    const out = Object.assign({}, layout || {});
    delete out.width;
    delete out.height;
    out.autosize = true;
    return out;
  }

  async function renderPlot(tile, group, key, seq) {
    const plotEl = tile.querySelector('.plot');
    tile.style.display = '';
    plotEl.innerHTML = '<p class="muted">Loading…</p>';
    const url = `/results/${jobId}/plot/${encodeURIComponent(group)}/${encodeURIComponent(key)}?seq=${seq}`;
    try {
      const resp = await fetch(url);
      if (resp.status === 404) {
        tile.style.display = 'none';
        return;
      }
      if (!resp.ok) {
        emptyTile(tile, `Failed to load (HTTP ${resp.status}).`);
        return;
      }
      const fig = await resp.json();
      if (typeof Plotly === 'undefined') {
        emptyTile(tile, 'Plotly library failed to load.');
        return;
      }
      plotEl.innerHTML = '';
      await Plotly.newPlot(
        plotEl, fig.data, normalizeLayout(fig.layout),
        {responsive: true, displaylogo: false},
      );
    } catch (e) {
      console.error('Plot render failed', key, group, e);
      emptyTile(tile, `Error: ${e.message || e}`);
    }
  }

  function renderStats(group) {
    const rows = (meta.stats && meta.stats[group]) || [];
    statsTable.innerHTML = rows.map(r =>
      `<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`
    ).join('');
  }

  function refreshSeqOptions(group) {
    const names = (meta.sequence_names && meta.sequence_names[group]) || [];
    seqSelect.innerHTML = names.map((n, i) =>
      `<option value="${i}">${escapeHtml(n)}</option>`
    ).join('');
    if (names.length > 1) {
      seqWrapper.style.display = '';
    } else {
      seqWrapper.style.display = 'none';
    }
  }

  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, c => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
    }[c]));
  }

  async function renderCombined() {
    const el = document.getElementById('plot-combined');
    if (!el) return;
    try {
      const resp = await fetch(`/results/${jobId}/plot/combined`);
      if (!resp.ok) return;
      const fig = await resp.json();
      Plotly.newPlot(el, fig.data, fig.layout, {responsive: true, displaylogo: false});
    } catch (e) { /* ignore */ }
  }

  async function refreshAll() {
    const group = groupSelect.value;
    const seq = parseInt(seqSelect.value || '0', 10);
    renderStats(group);
    await Promise.all(tiles.map(tile =>
      renderPlot(tile, group, tile.dataset.key, PER_REF.has(tile.dataset.key) ? seq : 0)
    ));
  }

  async function refreshPerRef() {
    const group = groupSelect.value;
    const seq = parseInt(seqSelect.value || '0', 10);
    await Promise.all(tiles
      .filter(t => PER_REF.has(t.dataset.key))
      .map(tile => renderPlot(tile, group, tile.dataset.key, seq))
    );
  }

  groupSelect.addEventListener('change', () => {
    refreshSeqOptions(groupSelect.value);
    refreshAll();
  });
  seqSelect.addEventListener('change', refreshPerRef);

  refreshSeqOptions(groupSelect.value);
  renderCombined();
  refreshAll();
})();
