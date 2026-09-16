"""Builds the option form, the settings of a run, and the command line from
--dump-options.

Each cmuts subcommand describes its arguments as JSON through the hidden
--dump-options flag. This module runs those dumps once, keeps the options a
user may tune, and turns the options of a job description or a settings file
into the settings of one run. The command line and the saved settings file
are two renderings of the same mapping, so they cannot disagree. Nothing
here names an individual option, so the form and the binary cannot drift
apart.
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
    """Returns the form widget for one option, or None if the form cannot
    render it.

    Flags become checkboxes, sets become groups of checkboxes, options with
    choices become selects, and numbers become number inputs. A free string
    names a file, which the form has no place for.
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


def none_choice(option: dict) -> str | None:
    """Returns the name that selects the empty set of a set option, or None if
    the option needs at least one choice."""
    return option.get("none_choice")


def default_set(option: dict) -> list[str]:
    """Returns the choices a set option starts with. A default of the none
    choice gives an empty list."""
    default = option["default"]
    if not default or default == none_choice(option):
        return []
    return default.split(",")


def default_value(option: dict):
    """Returns one option's default. A set option's default is a list of
    names."""
    if widget_of(option) == "set":
        return default_set(option)
    return option["default"]


def defaults_of(spec: dict) -> dict:
    """Returns the default of every option, keyed by name."""
    return {option["name"]: default_value(option) for option in spec["options"]}


# --- Options that depend on another option ---


def dependency(option: dict) -> dict | None:
    """Returns the option and choices one option depends on, or None if it
    always applies."""
    return option.get("applies_when")


def satisfied(option: dict, values: dict) -> bool:
    """Returns True if the option this one depends on is set to a choice it
    needs. An option that depends on nothing is always satisfied."""
    needs = dependency(option)

    return needs is None or values.get(needs["option"]) in needs["choices"]


def cmuts_version(specs: dict[str, dict]) -> str:
    return next(iter(specs.values()))["version"]


# --- The form ---


def annotate(sub: str, option: dict, defaults: dict) -> dict:
    """Returns a copy of one option with the extra fields the template needs.

    The label comes from the dump, or from the option name if the dump has
    none. A dependent option carries the field it depends on and the choices
    it needs, and starts hidden unless the defaults satisfy it. A required
    select keeps its null default, which the template shows as a placeholder
    the user must replace.
    """
    needs = dependency(option)
    out = dict(option)
    out["depends_field"] = field_name(sub, needs["option"]) if needs else None
    out["depends_choices"] = needs["choices"] if needs else []
    out["shown"] = satisfied(option, defaults)
    out["widget"] = widget_of(option)
    out["field"] = field_name(sub, option["name"])
    out["label"] = option.get("label") or option["name"].replace("-", " ").capitalize()
    out["step"] = "1" if option["type"] in INTEGER_TYPES else "any"
    out["choice_labels"] = option.get("choice_labels") or {}
    out["default_set"] = set(default_set(option)) if out["widget"] == "set" else set()
    out["accepts_none"] = none_choice(option) is not None
    return out


def checked_governor(spec: dict, option: dict) -> None:
    """Raises ValueError if an option depends on one the dump does not have,
    or on one with no choices. The form matches the value of that option
    against the choices, so it must have some."""
    needs = dependency(option)

    if needs is None:
        return

    governs = option_by_name(spec, needs["option"])

    if governs is None:
        raise ValueError(
            f"{option['name']}: depends on {needs['option']}, which is not an option")

    if not governs["choices"]:
        raise ValueError(
            f"{option['name']}: depends on {needs['option']}, which holds no choices")


def load_form_spec() -> tuple[dict[str, dict], list[dict], list[dict]]:
    """Runs every dump and builds the model the form renders.

    Returns three things: specs, the full dump of each subcommand; required,
    the exposed options a run must set, which the form shows up front; and
    sections, the rest as {group, options}, one per group, in dump order.
    """
    specs: dict[str, dict] = {}
    required: list[dict] = []
    sections: list[dict] = []
    for sub in SUBCOMMANDS:
        spec = dump_options(sub)
        specs[sub] = spec
        defaults = defaults_of(spec)
        groups: dict[str, list[dict]] = {}
        for option in exposed_options(spec):
            checked_governor(spec, option)
            if option["required"]:
                required.append(annotate(sub, option, defaults))
            else:
                groups.setdefault(option["group"], []).append(
                    annotate(sub, option, defaults))
        for group, options in groups.items():
            sections.append({"group": group, "options": options})
    return specs, required, sections


def option_place(sub: str, name: str) -> str:
    """Returns the name a message uses for one option."""
    return f"{sub}.{name}"


# --- The settings of one run ---


def option_setting(option: dict, sub: str, provided: dict):
    """Returns the value a run uses for one option: what the source gives, or
    the default. Raises ValueError if a required option has neither, which is
    where that requirement is enforced for every caller."""
    if option["name"] in provided:
        return provided[option["name"]]
    if option["required"]:
        raise ValueError(f"{option_place(sub, option['name'])}: a value is required")
    return default_value(option)


def run_settings(specs: dict[str, dict], provided: dict) -> dict:
    """Returns every exposed option's value for one run, keyed by subcommand
    and then by name. The mapping is complete, so a later cmuts with different
    defaults still replays the run."""
    return {
        sub: {
            option["name"]: option_setting(option, sub, provided.get(sub, {}))
            for option in exposed_options(spec)
        }
        for sub, spec in specs.items()
    }


def settings_document(specs: dict[str, dict], settings: dict) -> dict:
    """Returns the settings file for one run: the options it used and the
    version of cmuts that ran it."""
    return {VERSION_KEY: cmuts_version(specs), OPTIONS_KEY: settings}


# --- The command line ---


def matches_default(option: dict, value) -> bool:
    return value == default_value(option)


def set_text(option: dict, value: list[str]) -> str:
    """Returns the argument that one set option takes. The empty set becomes
    the name of the none choice."""
    return ",".join(value) if value else none_choice(option)


def value_text(option: dict, value) -> str:
    return set_text(option, value) if widget_of(option) == "set" else str(value)


def option_argument(option: dict, value) -> list[str]:
    """Returns the words one option adds to a command line, or nothing if its
    value is the default or a flag is off."""
    if value is None or matches_default(option, value):
        return []
    if widget_of(option) == "flag":
        return [f"--{option['name']}"] if value else []
    return [f"--{option['name']}", value_text(option, value)]


def option_args(spec: dict, settings: dict) -> list[str]:
    """Returns the extra arguments one subcommand gets: every exposed option
    set to something other than its default."""
    args: list[str] = []
    for option in exposed_options(spec):
        args.extend(option_argument(option, settings.get(option["name"])))
    return args


def all_option_args(specs: dict[str, dict], settings: dict) -> dict[str, list[str]]:
    return {sub: option_args(specs[sub], settings.get(sub, {})) for sub in specs}


# --- Reading a settings file ---


def option_by_name(spec: dict, name: str) -> dict | None:
    """Returns the named option of a dump, exposed or not, or None if the dump
    has no such option."""
    for option in spec["options"]:
        if option["name"] == name:
            return option
    return None


def option_named(spec: dict, name: str) -> dict | None:
    """Returns the named option if the form sets it, or None otherwise."""
    option = option_by_name(spec, name)

    return option if option is not None and is_exposed(option) else None


def choice_list(option: dict) -> str:
    return ", ".join(option["choices"])


def checked_flag(option: dict, value, where: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{where}: expected true or false")
    return value


def checked_set(option: dict, value, where: str) -> list[str]:
    """Returns the choices of one set option. A list that holds only the none
    choice also gives the empty set. Raises ValueError for the empty set if
    the option needs at least one choice."""
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
        raise ValueError(f"{where}: expected a list of names")
    if value == [none_choice(option)]:
        value = []
    if not value and none_choice(option) is None:
        raise ValueError(f"{where}: choose at least one of {choice_list(option)}")
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
    """Returns the options table of a settings file. Raises ValueError if the
    file or the table is not an object."""
    if not isinstance(document, dict):
        raise ValueError("the settings file must hold an object")
    table = document.get(OPTIONS_KEY, {})
    if not isinstance(table, dict):
        raise ValueError(f"{OPTIONS_KEY!r} must hold an object")
    return table


def checked_table(spec: dict, sub: str, table) -> dict:
    """Returns one subcommand's values from a settings file, checking every
    name, value and dependency against the dump."""
    if not isinstance(table, dict):
        raise ValueError(f"{sub}: must hold an object")
    values = {}
    for name, value in table.items():
        where = option_place(sub, name)
        option = option_named(spec, name)
        if option is None:
            raise ValueError(f"{where}: not an option the form sets")
        values[name] = checked_value(option, value, where)
    return checked_dependencies(spec, sub, values)


def dependency_text(needs: dict) -> str:
    """Returns the wording a message uses for one dependency."""
    return f"{needs['option']} set to {' or '.join(needs['choices'])}"


def checked_dependencies(spec: dict, sub: str, values: dict) -> dict:
    """Returns the values unchanged. Raises ValueError if one of them would
    reach the command line while the option it depends on is set to something
    else. A value equal to the default never reaches the command line, so it
    always passes, which is what lets a settings file list every option the
    form holds."""
    settings = defaults_of(spec) | values
    for name, value in values.items():
        option = option_named(spec, name)
        if matches_default(option, value) or satisfied(option, settings):
            continue
        where = option_place(sub, name)
        raise ValueError(
            f"{where}: applies only with {dependency_text(dependency(option))}")
    return values


def document_settings(specs: dict[str, dict], document) -> dict:
    """Returns the values a settings file gives, keyed by subcommand. Every
    name and value is checked, so a typo is refused rather than ignored."""
    settings = {}
    for sub, table in settings_options(document).items():
        if sub not in specs:
            raise ValueError(f"{sub}: not a step of the pipeline")
        settings[sub] = checked_table(specs[sub], sub, table)
    return settings


def settings_fields(specs: dict[str, dict], document) -> dict:
    """Returns the values a settings file sets, keyed by form field name."""
    return {
        field_name(sub, name): value
        for sub, values in document_settings(specs, document).items()
        for name, value in values.items()
    }
