from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint


def register_args(main_func):
    add_subcommand('prodigal', main_func, help='run genome annotation tool prodigal')
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
