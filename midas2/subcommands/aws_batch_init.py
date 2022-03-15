from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, command, backtick


def nvme_size_str():
    # When /mnt/nvme does not exist, do not crash (check=False); just return empty string.
    return backtick("""df -m /mnt/nvme | awk '{print $2}' | tail -1""", check=False)


def init_nvme(args):
    # TODO:  Generalize the magic numbers 838 and 1715518 (those are for AWS instance type r5.12xlarge).  # pylint: disable=fixme
    # https://github.com/czbiohub/MIDAS2.0/issues/17
    if nvme_size_str() != '1715518':
        # Raid, format, and mount the NVME drives attached to this instance.
        tsprint("Initializing instance NVME storage.")
        try:
            command("""set -o pipefail; lsblk | grep 838 | awk '{print "/dev/"$1}' | xargs -n 10 s3mi raid nvme""")
        except Exception as e:
            try:
                # Sometimes we've formatted it in a prior incarnation but the mountpoint can't exist in the container to tell us.
                # In those cases we can just try to mount it.
                command("""mount /dev/md0 /mnt/nvme""")
            except:
                raise e
        assert nvme_size_str() == '1715518', "Failed to initialize and mount instance NVME storage."
    else:
        tsprint("Instance NVME storage previously initialized.")
        if args.force:
            tsprint("Ignoring --force argument.  It is usually unnecessary to reinitialize AWS instance storage.")



def register_args(main_func):
    add_subcommand('aws_batch_init', main_func, help="init AWS Batch instance (never run outside AWS Batch)")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    uname = backtick("uname")
    assert uname == "Linux", f"Operating system {uname} is not Linux."
    init_nvme(args)
