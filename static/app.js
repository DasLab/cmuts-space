// Adds one condition row, fetched from the server with a fresh index so
// each row's fields stay distinct.
let conditionIndex = 0;

function addConditionRow() {
  conditionIndex += 1;
  fetch(`/condition-row?index=${conditionIndex}`)
    .then((r) => r.text())
    .then((html) => {
      document.getElementById("conditions").insertAdjacentHTML("beforeend", html);
    });
}

// Polls a running job's status and reloads the page when it settles. The
// running page marks its log element with the job id; polling starts from
// there, so no inline script depends on this file being loaded first.
function pollStatus(jobId) {
  fetch(`/results/${jobId}/status`)
    .then((r) => r.json())
    .then((s) => {
      document.getElementById("live-log").textContent = s.log || "";
      if (s.status === "done" || s.status === "error") {
        window.location.reload();
      } else {
        setTimeout(() => pollStatus(jobId), 1500);
      }
    })
    .catch(() => setTimeout(() => pollStatus(jobId), 3000));
}

const liveLog = document.getElementById("live-log");
if (liveLog && liveLog.dataset.job) {
  pollStatus(liveLog.dataset.job);
}

// Shows a message on a button, then restores its label.
function flashButton(button, message) {
  const original = button.textContent;
  button.textContent = message;
  setTimeout(() => { button.textContent = original; }, 1200);
}

function copyUrl(event) {
  const input = document.getElementById("result-url");
  input.select();
  navigator.clipboard.writeText(input.value).catch(() => document.execCommand("copy"));
  flashButton(event.currentTarget, "Copied");
}

// The report document, once its shell has loaded. The manifest node is the
// marker: the frame holds an empty placeholder document before that.
function loadedReportDocument() {
  const frame = document.getElementById("report-frame");
  const doc = frame && frame.contentDocument;
  return doc && doc.getElementById("cmuts-manifest") ? doc : null;
}

function saveBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

// A control's state lives in its properties, which serialized markup drops.
// These write one control's state into the attributes the markup keeps.
function writeCheckedAttribute(source, target) {
  if (source.checked) target.setAttribute("checked", "");
  else target.removeAttribute("checked");
}

function writeSelectedAttribute(source, target) {
  const options = target.querySelectorAll("option");
  for (let i = 0; i < options.length; i += 1) {
    if (i === source.selectedIndex) options[i].setAttribute("selected", "");
    else options[i].removeAttribute("selected");
  }
}

function writeControlState(source, target) {
  const box = source.type === "checkbox" || source.type === "radio";
  if (source.tagName === "SELECT") writeSelectedAttribute(source, target);
  else if (box) writeCheckedAttribute(source, target);
  else target.setAttribute("value", source.value);
}

// The report's markup, with every control's state written in, so the saved
// file opens showing the channels and references the page shows. A clone
// carries the attributes, which leaves the live page untouched. Both queries
// walk the same tree, so the two lists pair up element by element.
function reportSnapshot(doc) {
  const root = doc.documentElement.cloneNode(true);
  const sources = doc.querySelectorAll("input, select");
  const targets = root.querySelectorAll("input, select");
  for (let i = 0; i < sources.length; i += 1) writeControlState(sources[i], targets[i]);
  return `<!doctype html>\n${root.outerHTML}`;
}

// Saves the report as one standalone file. The report page inlines its own
// styles, plotly.js, and the figure data it has fetched, so its markup is a
// whole report: the views already opened stay interactive offline.
function downloadReport(event) {
  const button = event.currentTarget;
  const doc = loadedReportDocument();
  if (doc === null) {
    flashButton(button, "Report not ready");
    return;
  }
  const html = reportSnapshot(doc);
  saveBlob(new Blob([html], { type: "text/html" }), `${button.dataset.job}-report.html`);
}
