from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, backtick


def aws_batch_submit(args):
    pass


def register_args(main_func):
    add_subcommand('aws_batch_submit', main_func, help="submit aws batch command, keeping records in the s3 operaitons bucket")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    uname = backtick("uname")
    assert uname == "Linux", f"Operating system {uname} is not Linux."
    aws_batch_submit(args)
