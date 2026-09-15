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

// --- Submitting a job ---

// Every job goes to /run as a job description in the "job" field, with the
// uploaded files as parts that the description names by "uploads/<part>".
// The form builds its description from its fields. A bundled example is a
// description already, and it names its files under "examples/".

function showRunStatus(message, isError) {
  const box = document.getElementById("run-status");
  box.textContent = message;
  box.className = isError ? "error-msg" : "hint";
  box.hidden = message === "";
  if (message) box.scrollIntoView({ block: "nearest" });
}

// Disables the buttons that submit a job while one submission is in flight,
// so a second click does not submit a second job.
function setSubmitting(submitting) {
  document.querySelectorAll("[data-submits-job]").forEach((button) => {
    button.disabled = submitting;
  });
}

function showRunError(message) {
  showRunStatus(message, true);
  setSubmitting(false);
}

// Posts one job and opens its results page.
function postJob(body) {
  showRunStatus("Submitting the job. Large files can take a while to upload.", false);
  setSubmitting(true);
  fetch("/run", { method: "POST", body })
    .then((response) => response.json()
      .catch(() => ({ detail: `The server answered with status ${response.status}.` }))
      .then((data) => ({ response, data })))
    .then(({ response, data }) => {
      if (!response.ok) throw new Error(data.detail || "The server refused the job.");
      window.location.assign(data.url);
    })
    .catch((error) => showRunError(error.message));
}

// A page that the browser restores from its history keeps the state it had
// during the submission, so this clears that state.
window.addEventListener("pageshow", (event) => {
  if (!event.persisted || document.getElementById("run-status") === null) return;
  showRunStatus("", false);
  setSubmitting(false);
});

// Returns the request body for one job description and the file parts that
// it names, given as [part, file] pairs.
function jobBody(job, files) {
  const body = new FormData();
  body.append("job", JSON.stringify(job));
  files.forEach(([part, file]) => body.append(part, file));
  return body;
}

function runExample(name) {
  fetch(`/examples/${name}/job.json`)
    .then((response) => {
      if (!response.ok) throw new Error(`There is no example named ${name}.`);
      return response.json();
    })
    .then((job) => postJob(jobBody(job, [])))
    .catch((error) => showRunError(error.message));
}

// Adds one file as a part of the request, and returns the reference that
// the job description uses for it.
function attachFile(files, file) {
  const part = `file-${files.length}`;
  files.push([part, file]);
  return `uploads/${part}`;
}

// Returns one condition row as a condition of the job description, or null
// if the row has no treated reads.
function rowCondition(row, files) {
  const condition = { name: row.querySelector('input[type="text"]').value };
  row.querySelectorAll("input[data-role]").forEach((input) => {
    condition[input.dataset.role] = Array.from(input.files, (file) => attachFile(files, file));
  });
  return condition.treated.length ? condition : null;
}

// Returns the value that one option control sends, or undefined if the
// control leaves the option at its default. An unchecked box and an empty
// field leave it at its default.
function controlValue(control) {
  if (control.type === "checkbox") {
    if (!control.checked) return undefined;
    return control.closest(".set-field") ? control.value : true;
  }
  const raw = control.value.trim();
  if (raw === "") return undefined;
  return control.type === "number" ? Number(raw) : raw;
}

// Returns the options table of the job description. A field is named
// "opt.<subcommand>.<option>", and the boxes of a set option share one name.
function formOptions(form) {
  const table = {};
  form.querySelectorAll('[name^="opt."]').forEach((control) => {
    const value = controlValue(control);
    if (value === undefined) return;
    const [, sub, name] = control.name.split(".");
    table[sub] = table[sub] || {};
    if (control.closest(".set-field")) table[sub][name] = (table[sub][name] || []).concat(value);
    else table[sub][name] = value;
  });
  return table;
}

// Returns the job description that the form holds, and adds the files that
// it names to files.
function formJob(form, files) {
  const reference = attachFile(files, form.querySelector('input[name="fasta"]').files[0]);
  const conditions = Array.from(form.querySelectorAll(".condition-row"))
    .map((row) => rowCondition(row, files))
    .filter((condition) => condition !== null);
  return { reference, conditions, options: formOptions(form) };
}

function submitRunForm(event) {
  event.preventDefault();
  const files = [];
  const job = formJob(event.currentTarget, files);
  postJob(jobBody(job, files));
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
