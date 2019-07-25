# This example subcommand module is hidden by default.
#
# To expose it, just add the following line at the top of iggtools/__main__.py:
#
#     from iggtools.subcommands import example_subcommand # pylint: disable=unused-import
#
# Doing so will immediately expose "example_subcommand" as a subcommand with optional argument "-1".
# You can confirm this by running "python3 -m iggtools -h".
# Running "python3 -m iggtools example_subcommand" will execute the main() function below.
# Adding "-1" to that will set args.test_value to 1.
#
# To create your own new subcommand, simply clone this file, and import it.

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint


def register_args(main_func):
    subparser = add_subcommand('example_subcommand', main_func, help='run example_subcommand')
    # The "-1" argument below is only for this subcommand (and other subcommands
    # that expressly define it).  If you wish to expose it for every subcommand,
    # define it under argparser._add_shared_subcommand_args.
    subparser.add_argument('-1',
                           '--one',
                           action='store_const',
                           const=1,
                           dest='test_value',
                           help="set test_value to 1")
    return main_func


@register_args
def main(args):
    tsprint(f"Doing important work in subcommand {args.subcommand} with args {vars(args)}")
