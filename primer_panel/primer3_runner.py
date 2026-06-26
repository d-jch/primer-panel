"""Primer3 Boulder IO runner for Stage 2 primer design.

Stage 2 is a Primer3Plus-like design step.  It receives SEQUENCE_TEMPLATE
(design_template) and SEQUENCE_TARGET (extended_target relative coords) from
Stage 1, and returns primer pairs with Primer3 design metrics.

Stage 2 uses primer3plus-core to generate Boulder IO input matching
Primer3Plus web defaults.  Product size ranges come from Primer3Plus
settings, NOT from --target-size.  A coverage QC check verifies that
primer products cover SEQUENCE_TARGET (warning only, not a filter).
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from .config import PipelineConfig

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
# Primer result data class
# ------------------------------------------------------------------

@dataclass
class PrimerResult:
    """Result of Primer3 design for a single target.

    All coordinates are template-relative (design_template).
    Genomic product coordinates are computed by Stage 3 (in-silico PCR).
    """

    target_id: str
    primer_rank: int
    forward_primer: str
    reverse_primer: str
    forward_tm: float
    reverse_tm: float
    tm_diff: float
    forward_gc: float
    reverse_gc: float
    primer3_product_size: int        # Primer3-reported product size (template-relative)
    primer_pair_penalty: float
    primer_left_start: int
    primer_left_len: int
    primer_right_start: int
    primer_right_len: int
    primer3_status: str
    primer3_explain: str
    sequence_target_start_0based: int   # SEQUENCE_TARGET start (0-based)
    sequence_target_length: int         # SEQUENCE_TARGET length


# ------------------------------------------------------------------
# Primer3 availability check
# ------------------------------------------------------------------

def check_primer3_available(bin_path: str = "primer3_core") -> bool:
    """Check whether primer3_core is callable."""
    return shutil.which(bin_path) is not None


# ------------------------------------------------------------------
# Boulder IO input generation (via primer3plus-core adapter)
# ------------------------------------------------------------------


def _build_boulder_input(
    target_id: str,
    sequence: str,
    cfg: PipelineConfig,
    seq_target_start: int,
    seq_target_len: int,
) -> tuple[str, int]:
    """Build Primer3 Boulder IO input via primer3plus-core adapter.

    Uses Primer3Plus default settings from primer3plus-core package.
    User CLI overrides are passed only when explicitly specified.

    Returns:
        Tuple of (boulder_input_string, first_base_index).
    """
    from .primer3plus_core_adapter import build_primer3plus_input

    # Build overrides dict from cfg (only user-specified values)
    overrides = cfg.build_primer3_overrides()

    return build_primer3plus_input(
        sequence_id=target_id,
        sequence_template=sequence,
        sequence_target=f"{seq_target_start},{seq_target_len}",
        settings_file=cfg.primer3plus_settings_file,
        overrides=overrides,
    )


# ------------------------------------------------------------------
# Coverage validation
# ------------------------------------------------------------------

def check_target_coverage(
    primer_left_start: int,
    primer_left_len: int,
    primer_right_start: int,
    primer_right_len: int,
    seq_target_start_0based: int,
    seq_target_length: int,
    template_start: int,
) -> tuple[bool, str]:
    """QC check: verify that the primer product covers SEQUENCE_TARGET.

    This is a warning/bug-detector, NOT a filter.  Normal products should
    naturally cover the target because Primer3 was given SEQUENCE_TARGET.

    Coordinates are template-relative.  The product spans from
    PRIMER_LEFT_0 (left_start) to PRIMER_RIGHT_0 + 1 (right_start + 1).

    Returns (is_valid, detail_string).
    """
    product_start = primer_left_start
    product_end = primer_right_start + 1
    target_start = seq_target_start_0based
    target_end = seq_target_start_0based + seq_target_length

    errors: list[str] = []
    if product_start > target_start:
        errors.append(
            f"product_start({product_start}) > target_start({target_start})"
        )
    if product_end < target_end:
        errors.append(
            f"product_end({product_end}) < target_end({target_end})"
        )

    if errors:
        return False, "; ".join(errors)
    return True, ""


# ------------------------------------------------------------------
# Boulder IO output parsing
# ------------------------------------------------------------------

def _parse_boulder_output(raw: str) -> dict[str, str]:
    """Parse Primer3 Boulder IO output into a key-value dict."""
    result: dict[str, str] = {}
    for line in raw.strip().splitlines():
        if "=" in line:
            key, _, value = line.partition("=")
            result[key.strip()] = value.strip()
    return result


def _extract_primers(
    parsed: dict[str, str],
    target_id: str,
    num_return: int,
    seq_target_start: int,
    seq_target_len: int,
    first_base_index: int = 1,
) -> list[PrimerResult]:
    """Extract primer pairs from parsed Primer3 output.

    Primer3 returns coordinates in PRIMER_FIRST_BASE_INDEX base (usually 1).
    We convert all start positions to pipeline-internal 0-based:
        left_start_0based  = raw_left_start - first_base_index
        right_start_0based = raw_right_start - first_base_index

    PRIMER_RIGHT_i start is the rightmost base of the reverse complement on
    the template.  Product end (exclusive) = right_start_0based + 1.
    """
    results: list[PrimerResult] = []

    explain = parsed.get("PRIMER_PAIR_EXPLAIN", "")
    pair_count = int(parsed.get("PRIMER_PAIR_NUM_RETURNED", "0"))

    if pair_count == 0:
        # No primers found — return single record with status info
        results.append(PrimerResult(
            target_id=target_id,
            primer_rank=0,
            forward_primer="",
            reverse_primer="",
            forward_tm=0.0,
            reverse_tm=0.0,
            tm_diff=0.0,
            forward_gc=0.0,
            reverse_gc=0.0,
            primer3_product_size=0,
            primer_pair_penalty=0.0,
            primer_left_start=0,
            primer_left_len=0,
            primer_right_start=0,
            primer_right_len=0,
            primer3_status="no_primer",
            primer3_explain=explain,
            sequence_target_start_0based=seq_target_start,
            sequence_target_length=seq_target_len,
        ))
        return results

    for i in range(min(pair_count, num_return)):
        # Primer3 format: PRIMER_LEFT_{i}_SEQUENCE, PRIMER_RIGHT_{i}_TM, etc.
        left_seq = parsed.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
        right_seq = parsed.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
        left_tm = float(parsed.get(f"PRIMER_LEFT_{i}_TM", "0"))
        right_tm = float(parsed.get(f"PRIMER_RIGHT_{i}_TM", "0"))
        left_gc = float(parsed.get(f"PRIMER_LEFT_{i}_GC_PERCENT", "0"))
        right_gc = float(parsed.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", "0"))
        penalty = float(parsed.get(f"PRIMER_PAIR_{i}_PENALTY", "0"))
        product_size = int(parsed.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", "0"))

        # LEFT/RIGHT coords are "start,len" format in Primer3's base
        left_coords = parsed.get(f"PRIMER_LEFT_{i}", "0,0").split(",")
        right_coords = parsed.get(f"PRIMER_RIGHT_{i}", "0,0").split(",")
        raw_left_start = int(left_coords[0]) if len(left_coords) > 0 else 0
        left_len = int(left_coords[1]) if len(left_coords) > 1 else 0
        raw_right_start = int(right_coords[0]) if len(right_coords) > 0 else 0
        right_len = int(right_coords[1]) if len(right_coords) > 1 else 0

        # Convert Primer3 coordinates to pipeline-internal 0-based
        left_start = raw_left_start - first_base_index
        right_start = raw_right_start - first_base_index

        tm_diff = abs(left_tm - right_tm)

        results.append(PrimerResult(
            target_id=target_id,
            primer_rank=i + 1,
            forward_primer=left_seq,
            reverse_primer=right_seq,
            forward_tm=round(left_tm, 2),
            reverse_tm=round(right_tm, 2),
            tm_diff=round(tm_diff, 2),
            forward_gc=round(left_gc, 2),
            reverse_gc=round(right_gc, 2),
            primer3_product_size=product_size,
            primer_pair_penalty=round(penalty, 4),
            primer_left_start=left_start,
            primer_left_len=left_len,
            primer_right_start=right_start,
            primer_right_len=right_len,
            primer3_status="ok",
            primer3_explain="",
            sequence_target_start_0based=seq_target_start,
            sequence_target_length=seq_target_len,
        ))

    return results


# ------------------------------------------------------------------
# Public API: run Primer3 on one target
# ------------------------------------------------------------------

def run_primer3_for_target(
    target_id: str,
    sequence: str,
    cfg: PipelineConfig,
    seq_target_start: int,
    seq_target_len: int,
    return_input: bool = False,
) -> list[PrimerResult] | tuple[list[PrimerResult], str]:
    """Run Primer3 for a single target sequence (Primer3Plus-like).

    The *sequence* is the design_template (extended_target + primer_flank
    on each side).  Product size ranges use Primer3Plus defaults.

    Args:
        target_id: target identifier
        sequence: design_template sequence (SEQUENCE_TEMPLATE)
        cfg: pipeline config
        seq_target_start: extended_target start relative to template (0-based)
        seq_target_len: extended_target length
        return_input: if True, also return the Boulder IO input string

    Returns:
        list of PrimerResult, or (list, boulder_input) if return_input=True.

    Raises FileNotFoundError if primer3_core is not available.
    Raises RuntimeError if Primer3 exits with an error.
    """
    if not check_primer3_available(cfg.primer3_bin):
        raise FileNotFoundError(
            f"{cfg.primer3_bin} not found; install primer3 in the "
            f"primer_panel environment or pass --primer3-bin"
        )

    boulder_in, first_base_index = _build_boulder_input(
        target_id, sequence, cfg, seq_target_start, seq_target_len,
    )

    # Determine num_return for parsing (use override or Primer3Plus default 10)
    num_return = cfg.primer_num_return if cfg.primer_num_return is not None else 10

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".in", delete=False, prefix="p3_"
    ) as tmp:
        tmp.write(boulder_in)
        tmp_path = tmp.name

    try:
        proc = subprocess.run(
            [cfg.primer3_bin, tmp_path],
            capture_output=True,
            text=True,
            timeout=30,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Primer3 timed out for {target_id}")
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    if proc.returncode != 0:
        raise RuntimeError(
            f"Primer3 failed for {target_id} (exit {proc.returncode}): {proc.stderr}"
        )

    parsed = _parse_boulder_output(proc.stdout)
    results = _extract_primers(
        parsed, target_id, num_return, seq_target_start, seq_target_len,
        first_base_index=first_base_index,
    )

    if return_input:
        return results, boulder_in
    return results


# ------------------------------------------------------------------
# Sequence validation helpers
# ------------------------------------------------------------------

def is_all_n(sequence: str) -> bool:
    """Check if a sequence consists entirely of N characters."""
    return all(c in ("N", "n") for c in sequence if not c.isspace())


def has_n(sequence: str) -> bool:
    """Check if a sequence contains any N characters."""
    return any(c in ("N", "n") for c in sequence if not c.isspace())


def parse_fasta_sequences(fasta_path: Path) -> dict[str, str]:
    """Parse a FASTA file into {header_id: sequence} dict.

    The header_id is the first whitespace-delimited token after '>'.
    """
    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_seq_parts: list[str] = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq_parts)
                current_id = line[1:].split()[0]
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq_parts)

    return sequences
