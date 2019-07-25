from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, command, command_output, backtick


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



def register_args(main_func):
    add_subcommand('aws_batch_init', main_func, help="inialize AWS Batch instance")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    uname = backtick("uname")
    assert uname == "Linux", f"Operating system {uname} is not Linux."
    init_nvme(args)
