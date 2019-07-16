#!/usr/bin/env python3
#
# Import every subcommand you wish to expose.  Omitting a subcommand here will hide it.
#
# Pylint thinks these modules are "unused" but importing them has the side effect
# of registering their respective subcommands into the argument parser.
#
from iggtools.subcommands import example_a, example_b  # pylint: disable=unused-import


# ------------ no user-serviceable parts below this line -----------------------------
from iggtools.common import argparser


def main():
    # The one and only subcommand specified by the user will be invoked with its parsed args.
    args = argparser.parse_args()
    args.subcommand_main(args)


if __name__ == "__main__":
    main()
