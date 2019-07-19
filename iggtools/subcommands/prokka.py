from iggtools.common import argparser
from iggtools.common.utils import tsprint


def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand}.")
    tsprint(vars(args))
    assert args.subcommand == "prokka"


def register_argparser(singleton=[None]):  # pylint: disable=dangerous-default-value
    subparser = singleton[0]
    if not subparser:
        summary = 'run genome annotation tool prokka'
        subparser = argparser.get().subparsers.add_parser('prokka', description=summary, help=summary)
        subparser.set_defaults(subcommand_main=main)
        argparser.add_shared_subcommand_args(subparser)
        singleton[0] = subparser


register_argparser()
