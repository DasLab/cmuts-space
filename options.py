"""Derives the option form and the command-line arguments from --dump-options.

Each cmuts subcommand describes its arguments as JSON via the hidden
--dump-options flag. This module runs those dumps once, keeps the options a
user may tune, and turns submitted form values back into command-line
arguments. Nothing here names an individual option, so the form and the
binary cannot drift apart.
"""

from __future__ import annotations

import json
import shutil
import subprocess

# The subcommands whose options the form exposes, in pipeline order.
SUBCOMMANDS = ("align", "hmm", "sub", "div", "norm")

# The option groups the server owns. Paths, threads, and information flags
# never reach the form.
SERVER_GROUPS = {"Input and output", "Performance", "Information"}

NUMBER_TYPES = {"int", "size", "double"}
INTEGER_TYPES = {"int", "size"}


def cmuts_available() -> bool:
    return shutil.which("cmuts") is not None


def dump_options(sub: str) -> dict:
    """Runs one subcommand's --dump-options and parses the JSON it prints."""
    told = subprocess.run(
        ["cmuts", sub, "--dump-options"], capture_output=True, text=True,
    )
    if told.returncode != 0:
        raise RuntimeError(
            f"cmuts {sub} --dump-options failed: {told.stderr.strip() or 'refused'}"
        )
    return json.loads(told.stdout)


def widget_of(option: dict) -> str | None:
    """The form widget for one option, or None where it has no rendering.

    Flags become checkboxes, sets become checkbox groups, choice-restricted
    options become selects, and numbers become number inputs. A free string
    (a file path) has no place in the form.
    """
    if option["type"] == "flag":
        return "flag"
    if option["type"] == "set":
        return "set"
    if option["choices"]:
        return "select"
    if option["type"] in NUMBER_TYPES:
        return "number"
    return None


def is_exposed(option: dict) -> bool:
    return (
        not option["hidden"]
        and option["group"] not in SERVER_GROUPS
        and widget_of(option) is not None
    )


def field_name(sub: str, option: dict) -> str:
    return f"opt.{sub}.{option['name']}"


def annotate(sub: str, option: dict) -> dict:
    """One option extended with the fields the template reads directly. The
    dumped label names the field; the option name fills in where the dump
    carries none. A required select keeps its null default, which the
    template renders as a placeholder the user must replace."""
    out = dict(option)
    out["widget"] = widget_of(option)
    out["field"] = field_name(sub, option)
    out["label"] = option.get("label") or option["name"].replace("-", " ").capitalize()
    out["step"] = "1" if option["type"] in INTEGER_TYPES else "any"
    out["choice_labels"] = option.get("choice_labels") or {}
    default = option["default"]
    out["default_set"] = set(default.split(",")) if isinstance(default, str) else set()
    return out


def load_form_spec() -> tuple[dict[str, dict], list[dict], list[dict]]:
    """The dumps and the form model built from them.

    Returns (specs, required, sections): specs maps each subcommand to its
    full dump; required lists the exposed required options, which the form
    shows up front; sections lists the remaining exposed options as
    {group, options} in dump order, one section per option group.
    """
    specs: dict[str, dict] = {}
    required: list[dict] = []
    sections: list[dict] = []
    for sub in SUBCOMMANDS:
        spec = dump_options(sub)
        specs[sub] = spec
        groups: dict[str, list[dict]] = {}
        for option in spec["options"]:
            if not is_exposed(option):
                continue
            if option["required"]:
                required.append(annotate(sub, option))
            else:
                groups.setdefault(option["group"], []).append(annotate(sub, option))
        for group, options in groups.items():
            sections.append({"group": group, "options": options})
    return specs, required, sections


def parse_number(option: dict, raw: str):
    kind = int if option["type"] in INTEGER_TYPES else float
    try:
        return kind(raw)
    except ValueError as error:
        raise ValueError(f"--{option['name']}: {raw!r} is not a number") from error


def submitted_value(option: dict, sub: str, form) -> str | None:
    """One option's submitted value as command-line text, or None where the
    form leaves it at its default."""
    field = field_name(sub, option)
    widget = widget_of(option)

    if widget == "flag":
        # A flag can only be set on the command line, so an unchecked
        # default-true flag has no spelling and stays at its default.
        checked = form.get(field) is not None
        return "" if checked and not option["default"] else None

    if widget == "set":
        chosen = [c for c in option["choices"] if c in form.getlist(field)]
        if not chosen:
            return None
        value = ",".join(chosen)
        return value if value != option["default"] else None

    raw = (form.get(field) or "").strip()

    if widget == "select":
        if option["required"]:
            return raw or option["choices"][0]
        return raw if raw and raw != option["default"] else None

    if not raw:
        return None
    value = parse_number(option, raw)
    return str(value) if value != option["default"] else None


def option_args(sub: str, spec: dict, form) -> list[str]:
    """The extra arguments one subcommand receives from the form: each
    exposed option whose submitted value differs from its default."""
    args: list[str] = []
    for option in spec["options"]:
        if not is_exposed(option):
            continue
        value = submitted_value(option, sub, form)
        if value is None:
            continue
        args.append(f"--{option['name']}")
        if value != "":
            args.append(value)
    return args


def all_option_args(specs: dict[str, dict], form) -> dict[str, list[str]]:
    return {sub: option_args(sub, specs[sub], form) for sub in specs}
