#!/bin/bash
# Smoke-checks a running cmuts space: submits each bundled example dataset,
# waits for its pipeline to finish, and checks that its report answers.
#
#     scripts/smoke.sh BASE_URL [EXAMPLE...]
#
# EXAMPLE defaults to every bundled dataset. The same script covers the
# container built in CI and the live deployment; only BASE_URL differs.
#
# Author: Hamish M. Blair <hmblair@stanford.edu>

set -euo pipefail

WAIT_SECONDS=900
POLL_SECONDS=5

BASE=${1:?usage: smoke.sh BASE_URL [EXAMPLE...]}
shift
[ $# -gt 0 ] || set -- single-ref multi-ref

# Prints one field of a JSON body. Fails without a message where the body is
# not the expected JSON, so a call that does not answer prints no traceback.
json_field() {
    python3 -c 'import json, sys; print(json.load(sys.stdin)[sys.argv[1]])' "$1" 2>/dev/null
}

# Posts the job description of one example to /run, unchanged, and prints the
# URL of the job.
submit() {
    local url

    url=$(curl -sf --max-time 60 --retry 3 --retry-all-errors --retry-delay 2 \
        "$BASE/examples/$1/job.json" \
        | curl -sf --max-time 60 -F "job=<-" "$BASE/run" \
        | json_field url) || return 1
    echo "$BASE$url"
}

status_of() {
    curl -sf --max-time 30 "$1/status" | json_field status
}

report_answers() {
    curl -sf --max-time 60 --retry 3 --retry-all-errors --retry-delay 2 \
        -o /dev/null "$1/report/"
}

# Reports the failure and the job's pipeline log, then exits.
fail() {
    echo "smoke: $1" >&2
    curl -sf "$2/status" | json_field log >&2 || true
    exit 1
}

# Polls one job until it leaves the running state; prints the final status,
# or "timeout" where it never settles. The loop retries a poll that does not
# answer, so a server that stays unreachable reads as a timeout.
wait_until_done() {
    local deadline=$((SECONDS + WAIT_SECONDS))
    local state

    while [ "$SECONDS" -lt "$deadline" ]; do
        if state=$(status_of "$1"); then
            if [ "$state" != running ]; then
                echo "$state"
                return
            fi
        else
            echo "smoke: the status call did not answer; retrying" >&2
        fi
        sleep "$POLL_SECONDS"
    done
    echo timeout
}

check_example() {
    local job state

    job=$(submit "$1") || { echo "smoke: could not submit $1" >&2; exit 1; }
    echo "smoke: $1 -> $job"

    state=$(wait_until_done "$job")
    [ "$state" = done ] || fail "$1 finished as $state" "$job"
    report_answers "$job" || fail "$1 report does not answer" "$job"

    echo "smoke: $1 passed"
}

for example in "$@"; do
    check_example "$example"
done
