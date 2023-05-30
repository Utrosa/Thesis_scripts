####################################################################################################
# File:     cli.py
# Purpose:  Command line interface for the REPP package.
#
####################################################################################################
import click
from sing4me import __version__


@click.group()
@click.version_option(__version__, "--version", "-v", message="%(version)s")
def sing4me():
    pass
