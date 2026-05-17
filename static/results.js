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
    const options = [];
    if (names.length > 1) {
      // Special first option that maps to the multi-reference heatmap
      // view (the pre-saved profile.json) on the server.
      options.push('<option value="all">All references (heatmap)</option>');
    }
    names.forEach((n, i) => {
      options.push(`<option value="${i}">${escapeHtml(n)}</option>`);
    });
    seqSelect.innerHTML = options.join('');
    seqWrapper.style.display = names.length > 1 ? '' : 'none';
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

  // Tiles inside a closed <details> have zero size on first render, which
  // makes Plotly draw with width 0. Track which tiles are stale (need
  // rendering) and only render the ones currently visible. When a
  // <details> opens, render any stale tiles inside.
  const stale = new WeakSet();

  function isVisible(tile) {
    let el = tile.parentElement;
    while (el && el !== document.body) {
      if (el.tagName === 'DETAILS' && !el.open) return false;
      el = el.parentElement;
    }
    return true;
  }

  async function renderTileIfVisible(tile, group, seq) {
    if (!isVisible(tile)) {
      stale.add(tile);
      return;
    }
    stale.delete(tile);
    const key = tile.dataset.key;
    // Non-per-reference tiles always show the whole-dataset view.
    await renderPlot(tile, group, key, PER_REF.has(key) ? seq : '0');
  }

  async function refreshAll() {
    const group = groupSelect.value;
    const seq = seqSelect.value || '0';
    renderStats(group);
    // Mark every tile stale so closed-section tiles get re-rendered when
    // opened — important after group switches.
    tiles.forEach(t => stale.add(t));
    await Promise.all(tiles.map(t => renderTileIfVisible(t, group, seq)));
  }

  async function refreshPerRef() {
    const group = groupSelect.value;
    const seq = seqSelect.value || '0';
    await Promise.all(tiles
      .filter(t => PER_REF.has(t.dataset.key))
      .map(tile => renderTileIfVisible(tile, group, seq))
    );
  }

  function bindLazyDetails() {
    document.querySelectorAll('details').forEach(d => {
      d.addEventListener('toggle', async () => {
        if (!d.open) return;
        const group = groupSelect.value;
        const seq = seqSelect.value || '0';
        const inner = d.querySelectorAll('.plot-tile');
        await Promise.all(Array.from(inner).map(t => {
          if (stale.has(t)) return renderTileIfVisible(t, group, seq);
          // Already rendered — just nudge Plotly to refit if size changed.
          const plotEl = t.querySelector('.plot');
          if (typeof Plotly !== 'undefined' && plotEl && plotEl.firstChild) {
            try { Plotly.Plots.resize(plotEl); } catch (e) { /* ignore */ }
          }
        }));
      });
    });
  }

  groupSelect.addEventListener('change', () => {
    refreshSeqOptions(groupSelect.value);
    refreshAll();
  });
  seqSelect.addEventListener('change', refreshPerRef);

  bindLazyDetails();
  refreshSeqOptions(groupSelect.value);
  renderCombined();
  refreshAll();
})();
