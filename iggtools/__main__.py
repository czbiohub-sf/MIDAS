#!/usr/bin/env python3
#
# Import every subcommand you wish to expose.  Omitting a subcommand here will hide it.
#
# Pylint thinks these modules are "unused" but importing them has the side effect
# of registering their respective subcommands into the argument parser.
#
# Each subcommand module has its own main() function.  They are totally independent,
# with their own set of command line arguments and subcommand help text.
#
# The help text and broad description for the overall iggtools application is found in
# iggtools.common.argparser, and also there are defined any and all arguments that
# are shared across all subcommands.
#
from iggtools.subcommands import aws_batch_init, init, prokka  # pylint: disable=unused-import
from iggtools.common import argparser


def main():
    parser = argparser.get()
    args = parser.parse_args()
    # Invoke subcommand specified in args.
    subcommand_main = args.subcommand_main
    del args.subcommand_main  # Remove unserializable function pointer from args.
    return subcommand_main(args)


if __name__ == "__main__":
    main()
