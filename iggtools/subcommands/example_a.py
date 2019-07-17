# This example subcommand module is hidden by default.
#
# To expose it, just add the following line at the top of iggtools/__main__.py:
#
#     from iggtools.subcommands import example_a
#
# Doing so will immediately expose "example_a" as a subcommand with optional argument "-1".
# You can confirm this by running "python3 -m iggtools -h".
# Running "python3 -m iggtools example_a" will execute the main() function below.
# Adding "-1" to that will set args.test_value to 1.
#
# To create your own new subcommand, simply clone this file, and import it.


from iggtools.common import argparser
from iggtools.common.utils import tsprint


def main(args):
    tsprint(f"Doing important work in subcommand example_a.")
    tsprint(vars(args))
    assert args.subcommand == "example_a"


def register_argparser(singleton=[None]):  # pylint: disable=dangerous-default-value
    subparser = singleton[0]
    if not subparser:
        subparser = argparser.get().subparsers.add_parser('example_a', help='run subcommand example_a')
        subparser.add_argument('-1', '--one', action='store_const', const=1, dest='test_value', help="test_value := 1")
        subparser.set_defaults(subcommand_main=main)
        singleton[0] = subparser


register_argparser()
