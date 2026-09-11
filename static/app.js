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

function copyUrl(event) {
  const input = document.getElementById("result-url");
  input.select();
  navigator.clipboard.writeText(input.value).catch(() => document.execCommand("copy"));
  const btn = event.currentTarget;
  const orig = btn.textContent;
  btn.textContent = "Copied";
  setTimeout(() => { btn.textContent = orig; }, 1200);
}
