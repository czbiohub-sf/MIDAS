import json

from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, backtick, datecode, OutputStream
from iggtools.params.outputs import opsdir


def aws_batch_submit(args):
    name = args.batch_command.replace(" ", "__").replace('"', '_').replace("'", "_").replace('$', '_').replace('\\', '_').replace(":", "_").replace("-", "_")
    cmd = f"""aegea batch submit --name {name} --ecr-image {args.batch_ecr_image} --memory {args.batch_memory} --vcpus {args.batch_vcpus} --queue {args.batch_queue} --privileged --command="pip3 install 'git+https://github.com/czbiohub/iggtools.git@{args.batch_branch}' --upgrade ; iggtools --version ; aws s3 cp s3://microbiome-igg/2.0/README.TXT - ; iggtools aws_batch_init ; cd /mnt/nvme ; {args.batch_command} ; echo DONE" """
    aegea_output_json = backtick(cmd, quiet=False)
    ao = json.loads(aegea_output_json)
    job_id = ao['jobId']
    datestamp, timestamp = datecode().split("__")
    eventpath = f"{opsdir}/{datestamp}/{timestamp}__aws_batch_submit__{job_id}"
    with OutputStream(eventpath) as e:
        e.write(f"{timestamp}:{datestamp}:  Job ID:  {job_id}\n")
        e.write(f"{timestamp}:{datestamp}:  Job target:  {repr(args.batch_command)}\n")
        e.write(f"{timestamp}:{datestamp}:  Aegea command:  {repr(cmd)}\n")
    tsprint("You may watch the job with the command\n" + f"aegea batch watch {job_id}")


def register_args(main_func):
    subparser = add_subcommand('aws_batch_submit', main_func, help=f"submit command to run under AWS Batch, keeping records under {opsdir}")
    subparser.add_argument('--batch_command',
                           dest='batch_command',
                           required=False,
                           help="Command to run in AWS Batch container after installing iggtools")
    subparser.set_defaults(batch_command="iggtools --help")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    aws_batch_submit(args)
