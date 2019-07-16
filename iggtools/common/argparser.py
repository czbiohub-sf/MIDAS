from argparse import ArgumentParser
from iggtools import version


def get(singleton=[None]):  # pylint: disable=dangerous-default-value
    """Manage a single command line parser for all subcommands."""
    parser = singleton[0]
    if not parser:
        wiki = "https://github.com/czbiohub/iggtools/wiki"
        summary = f"Integrated Gut Genome Tools, Version {version}"
        parser = ArgumentParser(prog="iggtools", description=summary, epilog=f"For more information, see {wiki}")
        parser.add_argument('-v', '--version', action='version', version=summary)
        parser.subparsers = parser.add_subparsers(help="run specified subcommand", dest='subcommand', required=True)
        singleton[0] = parser
    return parser


def parse_args():
    return get().parse_args()
