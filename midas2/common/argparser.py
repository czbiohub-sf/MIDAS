from argparse import ArgumentParser, SUPPRESS, RawDescriptionHelpFormatter
from midas2 import version
from midas2.params import batch
from midas2.params.outputs import opsdir


# TODO:  Implement '--batch' in subcommands.  # pylint: disable=fixme
BATCH_FLAGS_IN_SUBCOMMANDS = False


def _new(_):
    """Create a single command line parser for all MIDAS subcommands."""
    wiki = "https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-2.0"
    summary = f"Metagenomic Intra-Species Diversity Analysis System 2.0 (MIDAS 2.0), Version {version}"
    parser = ArgumentParser(
        prog="midas2",
        description=summary,
        epilog=f"For more information, see {wiki}",
        formatter_class=RawDescriptionHelpFormatter
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


def _add_shared_subcommand_args(subparser, suppress_batch_args=False):
    """Arguments to be made available in every subcommand."""
    #
    # Adding the same flag separately to each subcommand parser makes it possible for
    # the end user to specify the flag after the subcommand name, like so
    #
    #     midas2 subcommand --flag      ----- THIS FEELS RIGHT
    #
    #  If we were to put the flag in the main argument parser above, the flag would be
    #  required to precede the subcommand, like so
    #
    #     midas2 --flag subcommand      ----- THIS FEELS UNNATURAL, DON'T DO IT
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
                           help="debug mode: skip cleanup on error, extra prints")
    subparser.set_defaults(debug=False)
    subparser.add_argument('--zzz_worker_mode',
                           action='store_const',
                           const=True,
                           dest='zzz_worker_mode',
                           help=SUPPRESS) # "worker mode is reserved for master-worker decomposition (do not use directly)"
    subparser.set_defaults(zzz_worker_mode=False)
    subparser.add_argument('--batch_branch',
                           dest='batch_branch',
                           required=False,
                           help=f"midas2 git branch for AWS batch, default '{batch.default_branch}'" if not suppress_batch_args else SUPPRESS)
    subparser.set_defaults(batch_branch=batch.default_branch)
    subparser.add_argument('--batch_memory',
                           dest='batch_memory',
                           required=False,
                           help=f"memory size for AWS batch container, default '{batch.default_memory}'" if not suppress_batch_args else SUPPRESS)
    subparser.set_defaults(batch_memory=batch.default_memory)
    subparser.add_argument('--batch_vcpus',
                           dest='batch_vcpus',
                           required=False,
                           help=f"vCPUs for AWS batch container, default '{batch.default_vcpus}'" if not suppress_batch_args else SUPPRESS)
    subparser.set_defaults(batch_vcpus=batch.default_vcpus)
    subparser.add_argument('--batch_queue',
                           dest='batch_queue',
                           required=False,
                           help=f"AWS batch queue, default '{batch.default_queue}'" if not suppress_batch_args else SUPPRESS)
    subparser.set_defaults(batch_queue=batch.default_queue)
    subparser.add_argument('--batch_ecr_image',
                           dest='batch_ecr_image',
                           required=False,
                           help=f"ecr image for AWS batch, default '{batch.default_ecr_image}'" if not suppress_batch_args else SUPPRESS)
    subparser.set_defaults(batch_ecr_image=batch.default_ecr_image)
    if BATCH_FLAGS_IN_SUBCOMMANDS:
        subparser.add_argument('-b',
                               '--batch',
                               action='store_const',
                               const=True,
                               dest='debug',
                               help=f"submit this command to AWS Batch with record keeping in {opsdir}; exposes additional options from the aws_batch_submit subcommand")
        subparser.set_defaults(batch=False)


# ----------------------------------------------------- #
#      No user-serviceable parts below this line.       #
# ----------------------------------------------------- #


def add_subcommand(subcommand_name, subcommand_main, *args, **kwargs):
    if 'help' in kwargs and 'description' not in kwargs:
        kwargs['description'] = kwargs['help']
    assert isinstance(singleton, ArgumentParser)
    subparser = singleton.subparsers.add_parser(subcommand_name, *args, **kwargs)
    subparser.set_defaults(subcommand_main=subcommand_main)
    _add_shared_subcommand_args(subparser, suppress_batch_args="aws_batch" not in subcommand_name)
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


if BATCH_FLAGS_IN_SUBCOMMANDS:
    _add_shared_subcommand_args(singleton)
