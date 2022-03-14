import json
import time

from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, backtick, datecode, OutputStream
from midas2.params.outputs import opsdir


def assert_have_aegea(min_version="3.2.1"):
    try:
        # Assert that aegea is installed and at least the minimum supported version.
        #
        # Tooling installed in the Dockerfile does not require these types of checks.
        #
        # A few of the midas2 admin subcommands (including this one) are supported to run
        # directly on a laptop or dev server, outside of docker, and (for extra operational
        # lightness and flexibility) even without being installed by a package manager --
        # with the downside of having to perform this check.  We should keep these checks
        # to a minimum.  If more creep up, we will require docker and erase these checks.
        #
        aegea, version = backtick("aegea --version | head -1").split()
        assert aegea == "aegea"
        vvv = tuple(int(v) for v in version.split("."))
        uuu = tuple(int(u) for u in min_version.split("."))
        assert vvv >= uuu, f"Aegea {version} is too old, please upgrade to {min_version} or above."
    except:
        tsprint("SUGGESTION:  Please 'pip3 install --upgrade aegea'")
        raise


def aws_batch_submit(args):
    """Submit given command to AWS Batch and log timestamped event under s3://operations/... folder in json format."""
    assert_have_aegea()
    # Replace anything that's not alphanumeric in batch_command with '_'
    name = str.join('', (c if c.isalnum() else '_' for c in args.batch_command))
    cmd = f"""aegea batch submit --name {name} --ecs-image {args.batch_ecr_image} --memory {args.batch_memory} --vcpus {args.batch_vcpus} --queue {args.batch_queue} --privileged --command="pip3 install 'git+https://github.com/czbiohub/MIDAS2.0.git@{args.batch_branch}' --upgrade ; midas2 --version ; aws s3 cp s3://microbiome-igg/2.0/README.TXT - ; midas2 aws_batch_init ; cd /mnt/nvme ; {args.batch_command} ; echo DONE" """
    tsprint(f"Submitting to AWS Batch queue {args.batch_queue}:  {args.batch_command}")
    aegea_output_json = backtick(cmd)
    ao = json.loads(aegea_output_json)
    job_id = ao['jobId']
    t_submit = int(time.time())
    datestamp, timestamp = datecode(t_submit).split("__")
    # timestamp is a string, and that's good, because JSON can lose resolution for large integers
    event = {
        "unix_timestamp": timestamp,
        "utc_date": datestamp,
        "type": "aws_batch_submit",
        "job_id": job_id,
        "job_target": args.batch_command,
        "aegea_command": cmd,
    }
    eventpath = f"{opsdir}/events/{datestamp}/{timestamp}__aws_batch_submit__{job_id}.json"
    with OutputStream(eventpath) as e:
        e.write(json.dumps(event))
    tsprint("You may watch the job with the command\n" + f"aegea batch watch {job_id}")


def register_args(main_func):
    subparser = add_subcommand('aws_batch_submit', main_func, help=f"submit command to run under AWS Batch, keeping records under {opsdir}")
    subparser.add_argument('--batch_command',
                           dest='batch_command',
                           required=False,
                           help="Command to run in AWS Batch container after installing midas2")
    subparser.set_defaults(batch_command="midas2 --help")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    aws_batch_submit(args)
