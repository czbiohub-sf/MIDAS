from iggtools.common import argparser
from iggtools.common.utils import tsprint, command, command_output


def nvme_size_str():
    # When /mnt/nvme does not exist, do not crash (check=False); just return empty string.
    return command_output("""df -m /mnt/nvme | awk '{print $2}' | tail -1""", check=False)


def init_nvme(args):
    # TODO:  Generalize the magic numbers 838 and 1715518 (those are for AWS instance type r5.12xlarge).
    if nvme_size_str() != '1715518':
        # Raid, format, and mount the NVME drives attached to this instance.
        tsprint("Initializing instance NVME storage.")
        command("""set -o pipefail; lsblk | grep 838 | awk '{print "/dev/"$1}' | xargs -n 10 s3mi raid nvme""")
        assert nvme_size_str() == '1715518', "Failed to initialize and mount instance NVME storage."
    else:
        tsprint("Instance NVME storage previously initialized.")
        if args.force:
            tsprint("Ignoring --force argument.  It is usually unnecessary to reinitialize AWS instance storage.")


def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand}.")
    tsprint(vars(args))
    assert args.subcommand == "aws_batch_init"
    init_nvme(args)


def register_argparser(singleton=[None]):  # pylint: disable=dangerous-default-value
    subparser = singleton[0]
    if not subparser:
        summary = 'inialize AWS Batch instance'
        subparser = argparser.get().subparsers.add_parser('aws_batch_init', description=summary, help=summary)
        subparser.set_defaults(subcommand_main=main)
        argparser.add_shared_subcommand_args(subparser)
        singleton[0] = subparser


register_argparser()
