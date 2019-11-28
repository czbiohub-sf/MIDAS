from argparse import ArgumentParser
from iggtools import version


def _new(_):
    """Create a single command line parser for all iggtools subcommands."""
    wiki = "https://github.com/czbiohub/iggtools/wiki"
    summary = f"Integrated Gut Genome Tools, Version {version}"
    parser = ArgumentParser(
        prog="iggtools",
        description=summary,
        epilog=f"For more information, see {wiki}"
    )
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=summary
    )
    parser.subparsers = parser.add_subparsers(
        help="specify subcommand",
        dest='subcommand',
        required=True
    )
    return parser


def _add_shared_subcommand_args(subparser):
    """Arguments to be made available in every subcommand."""
    #
    # Adding the same flag separately to each subcommand parser makes it possible for
    # the end user to specify the flag after the subcommand name, like so
    #
    #     iggtools subcommand --flag      ----- THIS FEELS RIGHT
    #
    #  If we were to put the flag in the main argument parser above, the flag would be
    #  required to precede the subcommand, like so
    #
    #     iggtools --flag subcommand      ----- THIS FEELS UNNATURAL, DON'T DO IT
    #
    #  In our experience, having flags before the subcommand feels like an inconsistency.
    #  We make sure ALL flags follow the subcommand.  This requires redefining every flag
    #  for each subcommand.
    #
    subparser.add_argument('-f',
                           '--force',
                           action='store_const',
                           const=True,
                           dest='force',
                           help="force rebuild of pre-existing outputs")
    subparser.set_defaults(force=False)
    subparser.add_argument('-g',
                           '--debug',
                           action='store_const',
                           const=True,
                           dest='debug',
                           help="debug mode, e.g. skip cleanup on error, extra prints, etc.")
    subparser.set_defaults(debug=False)


# ----------------------------------------------------- #
#      No user-serviceable parts below this line.       #
# ----------------------------------------------------- #


def add_subcommand(subcommand_name, subcommand_main, *args, **kwargs):
    if 'help' in kwargs and 'description' not in kwargs:
        kwargs['description'] = kwargs['help']
    assert isinstance(singleton, ArgumentParser)
    subparser = singleton.subparsers.add_parser(subcommand_name, *args, **kwargs)
    subparser.set_defaults(subcommand_main=subcommand_main)
    _add_shared_subcommand_args(subparser)
    return subparser


def parse_args():
    assert isinstance(singleton, ArgumentParser)
    subcommand_args = singleton.parse_args()
    subcommand_main = subcommand_args.subcommand_main
    del subcommand_args.subcommand_main  # Remove unserializable function pointer from args.
    return subcommand_main, subcommand_args


@_new
def singleton():
    assert False, "Through the magic of decorators, singleton becomes an instance of ArgumentParser."
