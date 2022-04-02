from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, find_files
from midas2.params.outputs import genomes

def init(args):
    if find_files(genomes("uhgg")):
        msg = f"TABLE_OF_CONTENTS exits for UHGG."
        tsprint(msg)
    if find_files(genomes("gtdb")):
        msg = f"TABLE_OF_CONTENTS exits for GTDB."
        tsprint(msg)


def register_args(main_func):
    add_subcommand('init', main_func, help=f"initialize target table-of-contents", epilog=init.__doc__)
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    init(args)
