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

// Every row takes a fresh index, so its fields are named apart from the rest.
let conditionIndex = 0;

// Adds one empty condition row and returns it once the server has answered.
function appendConditionRow() {
  conditionIndex += 1;

  return fetch(`/condition-row?index=${conditionIndex}`)
    .then((response) => response.text())
    .then((html) => {
      const holder = document.createElement("template");

      holder.innerHTML = html.trim();

      const row = holder.content.firstElementChild;

      document.getElementById("conditions").append(row);
      refreshRunButton();

      return row;
    });
}

function removeConditionRow(button) {
  button.closest(".condition-row").remove();
  refreshRunButton();
}

// --- Submitting a job ---

// Every job goes to /run as a job description in the "job" field, with the
// uploaded files as parts that the description names by "uploads/<part>".
// The form builds that description from its fields.

function showRunStatus(message, isError) {
  const box = document.getElementById("run-status");
  box.textContent = message;
  box.className = isError ? "error-msg" : "hint";
  box.hidden = message === "";
  if (message) box.scrollIntoView({ block: "nearest" });
}

// Disables the buttons that submit a job while one submission is in flight,
// so a second click does not submit a second job. A form that is not ready to
// run keeps its own button disabled afterwards.
function setSubmitting(submitting) {
  document.querySelectorAll("[data-submits-job]").forEach((button) => {
    button.disabled = submitting;
  });

  if (!submitting) refreshRunButton();
}

// Returns the name a message gives one control, taken from the label holding
// it.
function controlName(control) {
  const label = control.closest("label");
  const text = label && label.childNodes[0] ? label.childNodes[0].textContent : "";

  return text.replace("(required)", "").trim().toLowerCase() || control.name;
}

// Joins names as a sentence: "a and b", or "a, b and c".
function sentenceList(names) {
  if (names.length < 2) return names.join("");

  return `${names.slice(0, -1).join(", ")} and ${names[names.length - 1]}`;
}

// Names everything a run still needs, in the order the form asks for it: the
// reference, a name and treated reads for every condition on the page, and a
// value for every required option.
function missingForRun(form) {
  const fasta = form.querySelector('input[name="fasta"]');
  const rows = Array.from(form.querySelectorAll(".condition-row"));
  const missing = [];

  if (!fasta || !fasta.files.length) missing.push("a reference FASTA");

  if (!rows.length) missing.push("at least one condition");

  if (!rows.every(rowHasName)) missing.push("a name for every condition");

  if (!rows.every(rowHasTreated)) missing.push("treated reads in every condition");

  /* The conditions have checks of their own, which name a whole row rather
     than one of its fields. */
  form.querySelectorAll("[required]:invalid").forEach((control) => {
    if (control !== fasta && !control.closest(".condition-row")) {
      missing.push(controlName(control));
    }
  });

  return missing;
}

// Disables the run button until the form is ready. The tooltip of the wrapper
// names what is missing, and shows while the button is disabled.
function refreshRunButton() {
  const form = document.getElementById("run-form");

  if (!form) return;

  const button = form.querySelector('button[type="submit"]');
  const wrapper = button && button.closest(".submit-wrap");
  const missing = missingForRun(form);

  if (button) button.disabled = missing.length > 0;

  if (wrapper) {
    wrapper.title = missing.length ? `The run needs ${sentenceList(missing)}.` : "";
  }
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

// --- Loading a bundled example into the form ---

// Fetches one file the example names and returns it as a File, which is what a
// file input holds.
async function exampleFile(path) {
  const response = await fetch(`/${path}`);

  if (!response.ok) throw new Error(`${path} could not be read.`);

  const name = path.split("/").pop();

  return new File([await response.blob()], name);
}

// Puts files into one file input, which takes them only as a FileList.
function setFiles(input, files) {
  const holder = new DataTransfer();

  files.forEach((file) => holder.items.add(file));
  input.files = holder.files;
}

// Fills one input with the files at the paths given, if there are any.
async function fillFileInput(input, paths) {
  if (!paths || !paths.length) return;

  setFiles(input, await Promise.all(paths.map(exampleFile)));
}

// Builds one condition row from the example, with its name and the reads of
// each role it names.
async function fillConditionRow(condition) {
  const row = await appendConditionRow();

  row.querySelector('input[type="text"]').value = condition.name || "";

  for (const input of row.querySelectorAll("input[data-role]")) {
    await fillFileInput(input, condition[input.dataset.role]);
  }
}

// Returns the name the form gives one option, as options.field_name spells it.
function fieldName(sub, name) {
  return `opt.${sub}.${name}`;
}

// Writes the options of a job description into the form.
function applyJobOptions(options) {
  const fields = {};

  Object.entries(options || {}).forEach(([sub, table]) => {
    Object.entries(table).forEach(([name, value]) => {
      fields[fieldName(sub, name)] = value;
    });
  });

  applySettings(fields);
}

// Loads a bundled example into the form. The run is left for the user to
// start, so the example shows what a run is made of.
async function loadExample(name) {
  const form = document.getElementById("run-form");

  closeMenus();
  showRunStatus(`Loading the ${name} example…`, false);
  setSubmitting(true);

  try {
    const response = await fetch(`/examples/${name}/job.json`);

    if (!response.ok) throw new Error(`There is no example named ${name}.`);

    const job = await response.json();

    form.querySelectorAll(".condition-row").forEach((row) => row.remove());
    await fillFileInput(form.querySelector('input[name="fasta"]'), [job.reference]);

    for (const condition of job.conditions) {
      await fillConditionRow(condition);
    }

    applyJobOptions(job.options);
    setSubmitting(false);
    showRunStatus("The example is loaded. Change any option, then run it.", false);
  } catch (error) {
    showRunError(error.message);
  }
}

// Adds one file as a part of the request, and returns the reference that
// the job description uses for it.
function attachFile(files, file) {
  const part = `file-${files.length}`;
  files.push([part, file]);
  return `uploads/${part}`;
}

// Whether a row carries the name every condition needs.
function rowHasName(row) {
  return row.querySelector('input[type="text"]').value.trim().length > 0;
}

// Whether a row holds the treated reads every condition needs.
function rowHasTreated(row) {
  return row.querySelector('input[data-role="treated"]').files.length > 0;
}

// Returns the condition one row describes. Every row is sent, complete or not,
// so no condition is dropped without being reported.
function rowCondition(row, files) {
  const condition = { name: row.querySelector('input[type="text"]').value };
  row.querySelectorAll("input[data-role]").forEach((input) => {
    condition[input.dataset.role] = Array.from(input.files, (file) => attachFile(files, file));
  });
  return condition;
}

// Returns the value that one option control sends, or undefined if the
// control leaves the option at its default. An unchecked box and an empty
// field leave it at its default.
function controlValue(control) {
  if (control.type === "checkbox") return control.checked ? true : undefined;
  const raw = control.value.trim();
  if (raw === "") return undefined;
  return control.type === "number" ? Number(raw) : raw;
}

// Returns the names of the checked boxes of one set option. The list is
// empty when no box is checked.
function checkedChoices(field) {
  return Array.from(field.querySelectorAll("input:checked"), (box) => box.value);
}

// Adds one value to the options table. The field name has the form
// "opt.<subcommand>.<option>", which gives the place of the value in the table.
function addOption(table, fieldName, value) {
  const [, sub, name] = fieldName.split(".");
  table[sub] = table[sub] || {};
  table[sub][name] = value;
}

// Returns the options table of the job description. The table always holds
// the checked boxes of each set option, because an empty list is a valid setting.
function formOptions(form) {
  const table = {};
  form.querySelectorAll('[name^="opt."]').forEach((control) => {
    if (control.matches(":disabled") || control.closest(".set-field")) return;
    const value = controlValue(control);
    if (value !== undefined) addOption(table, control.name, value);
  });
  form.querySelectorAll(".set-field").forEach((field) => {
    const first = field.querySelector("input");
    if (first && !first.matches(":disabled")) addOption(table, first.name, checkedChoices(field));
  });
  return table;
}

// Returns the job description that the form holds, and adds the files that
// it names to files.
function formJob(form, files) {
  const reference = attachFile(files, form.querySelector('input[name="fasta"]').files[0]);
  const conditions = Array.from(form.querySelectorAll(".condition-row"))
    .map((row) => rowCondition(row, files));
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
  showDependentOptions();
  refreshRunButton();
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


// --- Options that depend on another option ---

// Whether the option a field depends on holds one of the choices it needs.
function dependencyMet(field) {
  const governing = document.querySelector(`[name="${field.dataset.dependsOn}"]`);
  if (!governing) return true;
  return field.dataset.dependsChoices.split(",").includes(governing.value);
}

// Shows every field whose dependency is met and hides the rest. A hidden field
// is disabled as well, so the browser leaves it out of its checks and the run
// never carries an option cmuts would refuse.
function showDependentOptions() {
  document.querySelectorAll("[data-depends-on]").forEach((field) => {
    const met = dependencyMet(field);
    field.hidden = !met;
    field.disabled = !met;
  });
}

// Checks a box again if the user unchecks the last box of a set option that
// needs at least one choice.
function keepOneChoice(event) {
  const field = event.target.closest(".set-field[data-needs-choice]");
  if (field && !field.querySelector("input:checked")) event.target.checked = true;
}

const runForm = document.getElementById("run-form");

if (runForm) {
  runForm.addEventListener("change", keepOneChoice);
  runForm.addEventListener("change", showDependentOptions);
  showDependentOptions();

  runForm.addEventListener("change", refreshRunButton);
  runForm.addEventListener("input", refreshRunButton);
  refreshRunButton();
}


// --- The download menu of one condition ---

// Closes every open menu, apart from the one the caller is using.
function closeMenus(inUse) {
  document.querySelectorAll("details.menu[open]").forEach((menu) => {
    if (menu !== inUse) menu.open = false;
  });
}

document.addEventListener("click", (event) => {
  const clicked = event.target.closest("details.menu");

  closeMenus(clicked);
});
