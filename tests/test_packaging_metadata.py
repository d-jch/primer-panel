"""Tests for package metadata configuration."""

from pathlib import Path
import subprocess
import sys
import tomllib


def test_pyproject_version_is_read_from_package_init():
    pyproject = tomllib.loads(Path("pyproject.toml").read_text(encoding="utf-8"))

    assert "version" not in pyproject["project"]
    assert "version" in pyproject["project"]["dynamic"]
    assert pyproject["tool"]["setuptools"]["dynamic"]["version"] == {
        "attr": "primer_panel.__version__",
    }


def test_project_builds_wheel(tmp_path):
    """Exercise setuptools metadata generation and wheel construction."""
    project_root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            str(project_root),
            "--no-deps",
            "--no-build-isolation",
            "--wheel-dir",
            str(tmp_path),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert list(tmp_path.glob("primer_panel-*.whl"))
