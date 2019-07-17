# This example subcommand module is hidden by default.
#
# To expose it, just add the following line at the top of iggtools/__main__.py:
#
#     from iggtools.subcommands import example_b
#
# Doing so will immediately expose "example_b" as a subcommand with optional argument "-2".
# You can confirm this by running "python3 -m iggtools -h".
# Running "python3 -m iggtools example_b" will execute the main() function below.
# Adding "-2" to that will set args.test_value to 2.
#
# To create your own new subcommand, simply clone this file, and import it.


from iggtools.common import argparser


def main(args):
    print("Doing important work in example subcommand B.")
    print(args)


def register_argparser(singleton=[None]):  # pylint: disable-dangerous-default
    subparser = singleton[0]
    if not subparser:
        subparser = argparser.get().subparsers.add_parser('example_b', help='run subcommand example_b')
        subparser.add_argument('-2', '--two', action='store_const', const=2, dest='test_value', help="test_value := 2")
        subparser.set_defaults(subcommand_main=main)
        singleton[0] = subparser


register_argparser()
