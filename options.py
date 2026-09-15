"""Derives the option form, the settings of a run, and the command-line
arguments from --dump-options.

Each cmuts subcommand describes its arguments as JSON via the hidden
--dump-options flag. This module runs those dumps once, keeps the options a
user may tune, and turns the options table of a job description or a
settings file into the settings of one run. The command line and the saved
settings file are two renderings of that one mapping, so they cannot
disagree. Nothing here names
an individual option, so the form and the binary cannot drift apart.
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

# The top-level keys of a settings file.
VERSION_KEY = "cmuts"
OPTIONS_KEY = "options"

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


def exposed_options(spec: dict) -> list[dict]:
    return [option for option in spec["options"] if is_exposed(option)]


def field_name(sub: str, name: str) -> str:
    return f"opt.{sub}.{name}"


def default_set(option: dict) -> list[str]:
    """The choices a set option holds by default."""
    default = option["default"]
    return default.split(",") if default else []


def default_value(option: dict):
    """One option's value where nothing sets it."""
    if widget_of(option) == "set":
        return default_set(option)
    return option["default"]


def cmuts_version(specs: dict[str, dict]) -> str:
    return next(iter(specs.values()))["version"]


# --- The form ---


def annotate(sub: str, option: dict) -> dict:
    """One option extended with the fields the template reads directly. The
    dumped label names the field; the option name fills in where the dump
    carries none. A required select keeps its null default, which the
    template renders as a placeholder the user must replace."""
    out = dict(option)
    out["widget"] = widget_of(option)
    out["field"] = field_name(sub, option["name"])
    out["label"] = option.get("label") or option["name"].replace("-", " ").capitalize()
    out["step"] = "1" if option["type"] in INTEGER_TYPES else "any"
    out["choice_labels"] = option.get("choice_labels") or {}
    out["default_set"] = set(default_set(option)) if out["widget"] == "set" else set()
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
        for option in exposed_options(spec):
            if option["required"]:
                required.append(annotate(sub, option))
            else:
                groups.setdefault(option["group"], []).append(annotate(sub, option))
        for group, options in groups.items():
            sections.append({"group": group, "options": options})
    return specs, required, sections


def option_place(sub: str, name: str) -> str:
    """How a message names one option."""
    return f"{sub}.{name}"


# --- The settings of one run ---


def option_setting(option: dict, sub: str, provided: dict):
    """One option's value for a run: what the source provides, or the
    option's default. A required option with neither raises, which is where
    the requirement is enforced for every caller."""
    if option["name"] in provided:
        return provided[option["name"]]
    if option["required"]:
        raise ValueError(f"{option_place(sub, option['name'])}: a value is required")
    return default_value(option)


def run_settings(specs: dict[str, dict], provided: dict) -> dict:
    """Every exposed option's value for one run, by subcommand. The mapping is
    complete, so a cmuts whose defaults have moved still replays the run."""
    return {
        sub: {
            option["name"]: option_setting(option, sub, provided.get(sub, {}))
            for option in exposed_options(spec)
        }
        for sub, spec in specs.items()
    }


def settings_document(specs: dict[str, dict], settings: dict) -> dict:
    """The settings file for one run: the options it used and the cmuts that
    ran it."""
    return {VERSION_KEY: cmuts_version(specs), OPTIONS_KEY: settings}


# --- The command line ---


def matches_default(option: dict, value) -> bool:
    return value == default_value(option)


def value_text(option: dict, value) -> str:
    return ",".join(value) if widget_of(option) == "set" else str(value)


def option_argument(option: dict, value) -> list[str]:
    """The command-line words one option contributes, empty where its value
    matches the default or a flag is off."""
    if value is None or matches_default(option, value):
        return []
    if widget_of(option) == "flag":
        return [f"--{option['name']}"] if value else []
    return [f"--{option['name']}", value_text(option, value)]


def option_args(spec: dict, settings: dict) -> list[str]:
    """The extra arguments one subcommand receives: each exposed option whose
    setting differs from its default."""
    args: list[str] = []
    for option in exposed_options(spec):
        args.extend(option_argument(option, settings.get(option["name"])))
    return args


def all_option_args(specs: dict[str, dict], settings: dict) -> dict[str, list[str]]:
    return {sub: option_args(specs[sub], settings.get(sub, {})) for sub in specs}


# --- Reading a settings file ---


def option_named(spec: dict, name: str) -> dict | None:
    for option in exposed_options(spec):
        if option["name"] == name:
            return option
    return None


def choice_list(option: dict) -> str:
    return ", ".join(option["choices"])


def checked_flag(option: dict, value, where: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{where}: expected true or false")
    return value


def checked_set(option: dict, value, where: str) -> list[str]:
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
        raise ValueError(f"{where}: expected a list of names")
    for name in value:
        if name not in option["choices"]:
            raise ValueError(f"{where}: {name!r} is not one of {choice_list(option)}")
    return value


def checked_choice(option: dict, value, where: str) -> str:
    if value not in option["choices"]:
        raise ValueError(f"{where}: {value!r} is not one of {choice_list(option)}")
    return value


def checked_number(option: dict, value, where: str):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{where}: expected a number")
    if option["type"] in INTEGER_TYPES and not isinstance(value, int):
        raise ValueError(f"{where}: expected a whole number")
    if option["minimum"] is not None and value < option["minimum"]:
        raise ValueError(f"{where}: {value} is below {option['minimum']}")
    if option["maximum"] is not None and value > option["maximum"]:
        raise ValueError(f"{where}: {value} is above {option['maximum']}")
    return value


CHECKS = {
    "flag": checked_flag,
    "set": checked_set,
    "select": checked_choice,
    "number": checked_number,
}


def checked_value(option: dict, value, where: str):
    return CHECKS[widget_of(option)](option, value, where)


def settings_options(document) -> dict:
    """The options table of a settings file."""
    if not isinstance(document, dict):
        raise ValueError("the settings file must hold an object")
    table = document.get(OPTIONS_KEY, {})
    if not isinstance(table, dict):
        raise ValueError(f"{OPTIONS_KEY!r} must hold an object")
    return table


def checked_table(spec: dict, sub: str, table) -> dict:
    """One subcommand's values from a settings file, with each name and value
    checked against the dump."""
    if not isinstance(table, dict):
        raise ValueError(f"{sub}: must hold an object")
    values = {}
    for name, value in table.items():
        where = option_place(sub, name)
        option = option_named(spec, name)
        if option is None:
            raise ValueError(f"{where}: not an option the form sets")
        values[name] = checked_value(option, value, where)
    return values


def document_settings(specs: dict[str, dict], document) -> dict:
    """The values a settings file provides, by subcommand. Every name and
    value is checked, so a typo cannot pass in silence."""
    settings = {}
    for sub, table in settings_options(document).items():
        if sub not in specs:
            raise ValueError(f"{sub}: not a step of the pipeline")
        settings[sub] = checked_table(specs[sub], sub, table)
    return settings


def settings_fields(specs: dict[str, dict], document) -> dict:
    """The form fields a settings file sets, keyed as the form names them."""
    return {
        field_name(sub, name): value
        for sub, values in document_settings(specs, document).items()
        for name, value in values.items()
    }
