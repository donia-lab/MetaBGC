"""Smoke tests for the metabgc command-line interface."""
import pytest
from click.testing import CliRunner

from metabgc import __version__
from metabgc.src.__main__ import cli

SUBCOMMANDS = [
    "build",
    "identify",
    "quantify",
    "cluster",
    "search",
    "analytics",
    "synthesize",
    "findtp",
]


def test_cli_help_lists_all_subcommands():
    result = CliRunner().invoke(cli, ["--help"])
    assert result.exit_code == 0
    for cmd in SUBCOMMANDS:
        assert cmd in result.output


def test_cli_version():
    result = CliRunner().invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


@pytest.mark.parametrize("cmd", SUBCOMMANDS)
def test_subcommand_help(cmd):
    result = CliRunner().invoke(cli, [cmd, "--help"])
    assert result.exit_code == 0
    assert "--help" in result.output
