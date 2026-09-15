// --- The job id in the embedding page's URL ---

// This app is served in a frame on its Hugging Face Space page. The URL of
// that page is the one the address bar shows and the one a reload asks for.
// The Space page hands its query string to the frame on load, so a job
// written there is the job a reload returns to.
const SPACE_PAGE_ORIGIN = "https://huggingface.co";

function writeJobIntoParentUrl(jobId) {
  if (window.parent === window) return;
  window.parent.postMessage(
    { queryString: jobId ? `?job=${jobId}` : "" },
    SPACE_PAGE_ORIGIN,
  );
}

// Every page runs this. A page with no job clears the query string, so
// returning to the form does not leave an old job in the URL.
const jobMarker = document.getElementById("job-marker");
writeJobIntoParentUrl(jobMarker ? jobMarker.dataset.job : "");

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

// --- Loading a settings file into the form ---

// Writes one value into one control, and reports whether the control moved.
// A set option shares its field name across several checkboxes, so a list
// checks each box the list names.
function applySetting(control, value) {
  const box = control.type === "checkbox";
  const before = box ? control.checked : control.value;
  if (box && Array.isArray(value)) control.checked = value.includes(control.value);
  else if (box) control.checked = Boolean(value);
  else control.value = String(value);
  return (box ? control.checked : control.value) !== before;
}

// Opens the sections holding a control, so a loaded setting is never hidden.
function revealControl(control) {
  let section = control.closest("details");
  while (section) {
    section.open = true;
    section = section.parentElement.closest("details");
  }
}

function applySettings(fields) {
  const form = document.getElementById("run-form");
  Object.keys(fields).forEach((name) => {
    form.querySelectorAll(`[name="${name}"]`).forEach((control) => {
      if (applySetting(control, fields[name])) revealControl(control);
    });
  });
}

function showSettingsError(message) {
  const box = document.getElementById("settings-error");
  box.textContent = message;
  box.hidden = message === "";
}

// Sends the chosen file to the server, which holds the option dumps, and
// writes the fields it answers with into the form. A refused file changes
// nothing, so the form never holds half of one.
function loadSettings(input) {
  const file = input.files[0];
  input.value = "";
  if (!file) return;
  const body = new FormData();
  body.append("settings", file);
  fetch("/settings", { method: "POST", body })
    .then((response) => response.json().then((data) => ({ response, data })))
    .then(({ response, data }) => {
      if (!response.ok) {
        showSettingsError(data.detail || "The settings file was refused.");
        return;
      }
      applySettings(data.fields);
      showSettingsError("");
    })
    .catch(() => showSettingsError("The settings file could not be read."));
}
