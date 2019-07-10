#!/usr/bin/env python3
import sys
assert sys.version_info >= (3, 7), "Python version >= 3.7 is required."

from argparse import ArgumentParser
from iggtools import version

def main():
    wiki = "https://github.com/czbiohub/iggtools/wiki"
    summary = f"Integrated Gut Genome Tools, Version {version}"
    parser = ArgumentParser(prog="iggtools", description=summary, epilog=f"For more information, see {wiki}")
    parser.add_argument('-v', '--version', action='version', version=summary)
    args = parser.parse_args()
    try:
        print(summary)
    except:
        parser.print_help()
        raise

if __name__ == "__main__":
    main()
