"""Dependency preflight checks for external bioinformatics tools.

Checks availability of primer3_core, isPcr, faToTwoBit, and validates
required input files.  Used by CLI stage gating and ``primer-panel doctor``.
"""

from __future__ import annotations

import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

# ── data structures ────────────────────────────────────────────────────────


@dataclass
class ToolCheck:
    """Result of checking one external tool."""

    name: str
    found: bool
    path: str | None = None
    install_hint: str = ""
    required: bool = True  # False = optional tool, won't affect all_ok


@dataclass
class FileCheck:
    """Result of checking one input file."""

    name: str
    path: Path | None
    exists: bool
    required: bool
    hint: str = ""


@dataclass
class DoctorReport:
    """Aggregated report from ``run_doctor``."""

    tools: list[ToolCheck]
    files: list[FileCheck]

    @property
    def all_ok(self) -> bool:
        tools_ok = all(t.found or not t.required for t in self.tools)
        files_ok = all(
            f.exists if f.path is not None else not f.required
            for f in self.files
        )
        return tools_ok and files_ok


# ── install hints ──────────────────────────────────────────────────────────

_PRIMER3_HINT = (
    "Primer3 is required for Stage 2 (primer design).\n"
    "  Install via micromamba/conda:\n"
    "    micromamba install -c bioconda primer3\n"
    "  Or pass --primer3-bin /path/to/primer3_core"
)

_ISPCR_HINT = (
    "UCSC isPcr is required for Stage 3 (specificity check).\n"
    "  Install via micromamba/conda:\n"
    "    micromamba install -c bioconda ispcr\n"
    "  Or pass --is-pcr-bin /path/to/isPcr"
)

_FATOTWOBIT_HINT = (
    "UCSC faToTwoBit is needed for --prepare-ispcr-db.\n"
    "  Install via micromamba/conda:\n"
    "    micromamba install -c bioconda ucsc-fatotwobit"
)


# ── tool checks ────────────────────────────────────────────────────────────


def check_tool(name: str, bin_path: str, install_hint: str, *, required: bool = True) -> ToolCheck:
    """Check whether *bin_path* is resolvable via ``shutil.which``."""
    resolved = shutil.which(bin_path)
    return ToolCheck(
        name=name,
        found=resolved is not None,
        path=resolved,
        install_hint=install_hint,
        required=required,
    )


def check_primer3(bin_path: str = "primer3_core") -> ToolCheck:
    return check_tool("primer3_core", bin_path, _PRIMER3_HINT, required=True)


def check_ispcr(bin_path: str = "isPcr") -> ToolCheck:
    return check_tool("isPcr", bin_path, _ISPCR_HINT, required=True)


def check_fatotwobit(bin_path: str = "faToTwoBit") -> ToolCheck:
    return check_tool("faToTwoBit", bin_path, _FATOTWOBIT_HINT, required=False)


# ── file checks ────────────────────────────────────────────────────────────


def check_file(
    name: str,
    path: Path | None,
    *,
    required: bool,
    hint: str = "",
) -> FileCheck:
    exists = path is not None and path.is_file()
    return FileCheck(name=name, path=path, exists=exists, required=required, hint=hint)


# ── preflight gates (called before each stage) ─────────────────────────────


def _validate_fasta_path(genome_fasta: Path | None, context: str) -> None:
    """Validate genome FASTA path exists and is a regular file.  Exits on failure."""
    if genome_fasta is None:
        _abort(
            f"--genome-fasta is required for {context}.\n"
            "  Provide a bgzipped+indexed hg38 FASTA:\n"
            "    primer-panel --genes ... --genome-fasta /path/to/hg38.fa"
        )
    if not genome_fasta.exists():
        _abort(
            f"--genome-fasta path not found: {genome_fasta}\n"
            f"  The FASTA file is required for {context}.\n"
            "  Check that the path is correct and the file exists."
        )
    if not genome_fasta.is_file():
        _abort(
            f"--genome-fasta is not a regular file: {genome_fasta}\n"
            "  Provide a bgzipped+indexed hg38 FASTA file."
        )


def preflight_stage2(cfg) -> None:
    """Validate prerequisites for Stage 2.  Exits on failure."""
    tc = check_primer3(cfg.primer3_bin)
    if not tc.found:
        _abort(tc.install_hint)


def preflight_stage3(cfg) -> None:
    """Validate prerequisites for Stage 3.  Exits on failure."""
    tc = check_ispcr(cfg.is_pcr_bin)
    if not tc.found:
        _abort(tc.install_hint)

    _validate_fasta_path(cfg.genome_fasta, "Stage 3 (specificity check)")


def preflight_prepare_ispcr_db(cfg) -> None:
    """Validate faToTwoBit is available for --prepare-ispcr-db."""
    tc = check_fatotwobit()
    if not tc.found:
        _abort(tc.install_hint)


def preflight_genome_fasta(genome_fasta: Path | None, stage: str | None) -> None:
    """Validate genome-fasta for the resolved stage.

    Checks:
    - stage=all (explicit or default): genome-fasta required and must exist
    - stage=design: genome-fasta required and must exist (Primer3 needs real sequences)
    - stage=specificity: genome-fasta required and must exist
    - stage=targets: no check (genome-fasta optional)

    Default (no --stage) is equivalent to --stage all.
    """
    # stage=targets does not require genome-fasta
    if stage == "targets":
        return
    # All other stages (design, specificity, all, or None/default) require it
    context = f"stage={stage}" if stage else "stage=all (the default pipeline)"
    _validate_fasta_path(genome_fasta, context)


# ── doctor (full diagnostic report) ───────────────────────────────────────


def run_doctor(
    *,
    genome_fasta: Path | None = None,
    common_dbsnp_bed: Path | None = None,
    primer3_bin: str = "primer3_core",
    ispcr_bin: str = "isPcr",
) -> DoctorReport:
    """Run all checks and return a structured report (never exits)."""
    tools = [
        check_primer3(primer3_bin),
        check_ispcr(ispcr_bin),
        check_fatotwobit(),
    ]
    files = [
        check_file(
            "genome FASTA",
            genome_fasta,
            required=False,
            hint="Required for Stage 2 real sequences and Stage 3 specificity.",
        ),
        check_file(
            "common dbSNP BED",
            common_dbsnp_bed,
            required=False,
            hint="Optional. Annotates primers with common variant overlap.",
        ),
    ]
    return DoctorReport(tools=tools, files=files)


def print_doctor_report(report: DoctorReport) -> None:
    """Pretty-print a DoctorReport to stderr."""
    print("=== primer-panel doctor ===", file=sys.stderr)
    print("", file=sys.stderr)

    print("External tools:", file=sys.stderr)
    for tc in report.tools:
        status = "OK" if tc.found else "MISSING"
        req_tag = "" if tc.required else " (optional)"
        loc = f" ({tc.path})" if tc.path else ""
        print(f"  [{status}] {tc.name}{req_tag}{loc}", file=sys.stderr)
        if not tc.found and tc.install_hint:
            for line in tc.install_hint.splitlines():
                print(f"         {line}", file=sys.stderr)

    print("", file=sys.stderr)
    print("Input files:", file=sys.stderr)
    for fc in report.files:
        if fc.path is None:
            status = "not provided"
        elif fc.exists:
            status = "OK"
        else:
            status = "NOT FOUND"
        tag = " (required)" if fc.required else " (optional)"
        print(f"  [{status}] {fc.name}{tag}", file=sys.stderr)
        if not fc.exists and fc.hint:
            for line in fc.hint.splitlines():
                print(f"         {line}", file=sys.stderr)

    print("", file=sys.stderr)
    if report.all_ok:
        print("All checks passed.", file=sys.stderr)
    else:
        print("Some checks failed — see above for details.", file=sys.stderr)


# ── internal ───────────────────────────────────────────────────────────────


def _abort(message: str) -> None:
    """Print an error message and exit with code 1."""
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)
