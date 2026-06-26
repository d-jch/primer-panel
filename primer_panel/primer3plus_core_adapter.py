"""Adapter for primer3plus-core package.

Bridges primer3plus-core's settings/boulder modules with our pipeline,
producing Primer3 Boulder IO input that matches Primer3Plus web defaults.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from primer3plus_core import boulder as p3p_boulder
from primer3plus_core import settings as p3p_settings

# Tags to ignore (Primer3Plus internal / UI-only / handled separately)
_IGNORE_TAGS = {
    "PRIMER_EXPLAIN_FLAG",
    "PRIMER_THERMODYNAMIC_PARAMETERS_PATH",
    "PRIMER_MIN_THREE_PRIME_DISTANCE",
    "P3P_SERVER_SETTINGS_FILE",
    "PRIMER_MASK_3P_DIRECTION",
    "PRIMER_MASK_5P_DIRECTION",
    "PRIMER_MASK_FAILURE_RATE",
    "PRIMER_MASK_KMERLIST_PATH",
    "PRIMER_MASK_KMERLIST_PREFIX",
    "PRIMER_MASK_TEMPLATE",
    "PRIMER_WT_MASK_FAILURE_RATE",
}

# P3_FILE_* tags (Primer3 file format metadata, not runtime tags)
_P3_FILE_PREFIXES = ("P3_FILE_", "P3_COMMENT", "P3_FILE_FLAG")

# P3P_* interface-only tags (always removed)
_P3P_PREFIX = "P3P_"

# Tags that should be removed if empty/placeholder
_REMOVE_IF_EMPTY = {
    "PRIMER_MISPRIMING_LIBRARY",
    "PRIMER_INTERNAL_MISHYB_LIBRARY",
    "PRIMER_MUST_MATCH_FIVE_PRIME",
    "PRIMER_MUST_MATCH_THREE_PRIME",
    "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME",
    "PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME",
}


def load_settings_dict(settings_file: Path | None = None) -> dict[str, str]:
    """Load Primer3Plus settings as a flat {tag: value} dict.

    Args:
        settings_file: Path to a settings JSON or .txt file.
            If None, uses primer3plus-core bundled default_settings.json.

    Returns:
        Flat dict of Primer3 tag → string value (only "def" entries).
    """
    if settings_file is None:
        raw = p3p_settings.load_default_settings()
    else:
        raw = settings_file.read_text(encoding="utf-8")

    # Try JSON first
    try:
        data = json.loads(raw)
        defs = data.get("def", data)
        return {k: v[0] if isinstance(v, list) else str(v) for k, v in defs.items()}
    except (json.JSONDecodeError, KeyError, TypeError):
        pass

    # Fall back to Boulder-IO text format (key=value per line)
    result: dict[str, str] = {}
    for line in raw.splitlines():
        line = line.strip()
        # Skip headers, comments, empty lines
        if not line or line.startswith("#") or line.startswith("Primer3 File"):
            continue
        if "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        val = val.strip()
        # Skip P3_FILE_* / P3_* non-Primer3 tags
        if key.startswith("P3_FILE_") or key in ("P3_COMMENT", "P3_FILE_FLAG", "P3_FILE_ID", "P3_FILE_TYPE"):
            continue
        if key:
            result[key] = val
    return result


def get_first_base_index(settings_file: Path | None = None,
                         overrides: dict[str, str] | None = None) -> int:
    """Get PRIMER_FIRST_BASE_INDEX from settings/overrides.

    Primer3Plus default is 1.
    """
    if overrides and "PRIMER_FIRST_BASE_INDEX" in overrides:
        return int(overrides["PRIMER_FIRST_BASE_INDEX"])
    settings = load_settings_dict(settings_file)
    return int(settings.get("PRIMER_FIRST_BASE_INDEX", "1"))


def build_primer3plus_input(
    sequence_id: str,
    sequence_template: str,
    sequence_target: str,
    settings_file: Path | None = None,
    overrides: dict[str, str] | None = None,
    first_base_index: int | None = None,
) -> tuple[str, int]:
    """Build Primer3 Boulder IO input matching Primer3Plus defaults.

    Args:
        sequence_id: SEQUENCE_ID value.
        sequence_template: Raw template sequence (cleaned to ACGTN).
        sequence_target: SEQUENCE_TARGET value in 0-based, e.g. "500,2700".
            Converted to Primer3's PRIMER_FIRST_BASE_INDEX for output.
        settings_file: Optional path to settings JSON/txt.
        overrides: User-specified tag overrides (from CLI args).
        first_base_index: Override for PRIMER_FIRST_BASE_INDEX (default: from settings, usually 1).

    Returns:
        Tuple of (boulder_input_string, first_base_index_used).
    """
    settings = load_settings_dict(settings_file)

    # Determine first_base_index
    if first_base_index is None:
        first_base_index = get_first_base_index(settings_file, overrides)

    # Clean template
    sequence_template = re.sub(r"[^ACGTNacgtn]", "", sequence_template).upper()

    # Build working dict
    tags: dict[str, str] = {}

    # Start from settings (skip P3P_*, P3_FILE_*, and ignored tags)
    for key, val in settings.items():
        if key.startswith(_P3P_PREFIX):
            continue
        if any(key.startswith(p) for p in _P3_FILE_PREFIXES):
            continue
        if key in _IGNORE_TAGS:
            continue
        if key in _REMOVE_IF_EMPTY and (not val or val.upper() == "NONE"):
            continue
        tags[key] = val

    # Remove empty SEQUENCE_* tags (except the ones we set below)
    for key in list(tags.keys()):
        if key.startswith("SEQUENCE_") and key not in (
            "SEQUENCE_TEMPLATE", "SEQUENCE_TARGET", "SEQUENCE_ID",
        ):
            if not tags[key] or tags[key] in ("", "-1000000", "-2000000"):
                del tags[key]

    # Convert 0-based sequence_target to Primer3's first_base_index
    # Pipeline internal: 0-based → "start,length"
    # Primer3 input: (start + first_base_index),length
    target_parts = sequence_target.split(",")
    if len(target_parts) == 2:
        target_start_0based = int(target_parts[0])
        target_len = int(target_parts[1])
        target_start_for_primer3 = target_start_0based + first_base_index
        sequence_target_for_primer3 = f"{target_start_for_primer3},{target_len}"
    else:
        sequence_target_for_primer3 = sequence_target

    # Set our sequence values
    tags["SEQUENCE_ID"] = sequence_id
    tags["SEQUENCE_TEMPLATE"] = sequence_template
    tags["SEQUENCE_TARGET"] = sequence_target_for_primer3

    # Always set EXPLAIN_FLAG=1 (Primer3Plus behavior)
    tags["PRIMER_EXPLAIN_FLAG"] = "1"

    # Apply user overrides
    if overrides:
        for key, val in overrides.items():
            if val is not None:
                tags[key] = str(val)

    # Build output with stable ordering:
    # SEQUENCE_ID, SEQUENCE_TEMPLATE, SEQUENCE_TARGET first, then sorted rest
    sequence_keys = ["SEQUENCE_ID", "SEQUENCE_TEMPLATE", "SEQUENCE_TARGET"]
    other_keys = sorted(k for k in tags if k not in sequence_keys)

    lines: list[str] = []
    for key in sequence_keys:
        if key in tags:
            lines.append(f"{key}={tags[key]}")
    for key in other_keys:
        lines.append(f"{key}={tags[key]}")
    lines.append("=")

    raw_input = "\n".join(lines) + "\n"

    # Let primer3plus-core do its preparation (strip thermo path, inject mispriming libs)
    prepared = p3p_boulder.prepare_input(raw_input)
    return prepared, first_base_index
