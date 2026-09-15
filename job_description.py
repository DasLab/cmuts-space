"""Reads the description of one job.

A job description is a JSON object. It names the reference FASTA and the
read files of each condition by role, and it holds the options of the run in
the format of a settings file:

    {
      "reference": "examples/single-ref/ref.fasta",
      "conditions": [
        {"name": "2A3", "treated": ["uploads/file-0"], "untreated": ["uploads/file-1"]}
      ],
      "options": {"align": {"preset": "map-ont"}}
    }

Each file is a reference that starts with its source. "uploads/<part>" names
a file part of the request, and "examples/<path>" names a bundled example
file. This module checks the shape of a description. The app resolves the
references.
"""

from __future__ import annotations

from dataclasses import dataclass

from pipeline import CONDITION_ROLES, MAX_CONDITIONS, ConditionInput

# A role holds a read file and, for paired-end reads, its mate.
MAX_FILES_PER_ROLE = 2


@dataclass
class JobDescription:
    """The reference and conditions of one job, with each file given as a
    reference that the app has not yet resolved."""

    reference: str
    conditions: list[ConditionInput]


def checked_reference(value) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError("reference: a reference FASTA is required")
    return value


def checked_files(value, where: str) -> list[str]:
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
        raise ValueError(f"{where}: expected a list of file references")
    if len(value) > MAX_FILES_PER_ROLE:
        raise ValueError(f"{where}: at most {MAX_FILES_PER_ROLE} files "
                         "(a read file and its mate)")
    return value


def checked_condition(index: int, condition) -> ConditionInput:
    """Returns one condition. An empty name becomes the condition's place in
    the list."""
    where = f"condition {index + 1}"
    if not isinstance(condition, dict):
        raise ValueError(f"{where}: expected an object")
    unknown = sorted(set(condition) - {"name", *CONDITION_ROLES})
    if unknown:
        raise ValueError(f"{where}: unknown fields {', '.join(unknown)}")
    name = condition.get("name", "")
    if not isinstance(name, str):
        raise ValueError(f"{where}.name: expected a string")
    files = {role: checked_files(condition.get(role, []), f"{where}.{role}")
             for role in CONDITION_ROLES}
    if not files["treated"]:
        raise ValueError(f"{where}: treated reads are required")
    return ConditionInput(name=name.strip() or where, **files)


def checked_conditions(value) -> list[ConditionInput]:
    if not isinstance(value, list) or not value:
        raise ValueError("conditions: at least one condition is required")
    if len(value) > MAX_CONDITIONS:
        raise ValueError(f"conditions: at most {MAX_CONDITIONS} conditions per run")
    return [checked_condition(i, c) for i, c in enumerate(value)]


def read_job_description(document) -> JobDescription:
    """Returns the reference and conditions of a job description. Raises
    ValueError if the description does not have the expected shape."""
    if not isinstance(document, dict):
        raise ValueError("the job description must be an object")
    return JobDescription(
        reference=checked_reference(document.get("reference")),
        conditions=checked_conditions(document.get("conditions")),
    )
