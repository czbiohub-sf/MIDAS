#!/usr/bin/env python3

import sys
import time
import subprocess
import multiprocessing


# Thread-safe and timestamped prints.
tslock = multiprocessing.RLock()


def timestamp(t):
    # We do not use "{:.3f}".format(time.time()) because its result may be
    # up to 0.5 ms in the future due to rounding.  We want flooring here.
    s = str(int(t * 10))
    return s[:-1] + "." + s[-1:]


def tsfmt(msg):
    ts = timestamp(time.time()) + ":  "
    msg = ts + str(msg).replace("\n", "\n" + ts)
    return msg


def tsout(msg):
    with tslock:
        sys.stdout.write(str(msg))
        sys.stdout.write("\n")


def tserr(msg):
    with tslock:
        sys.stderr.write(str(msg))
        sys.stderr.write("\n")


def tsprint(msg):
    tserr(tsfmt(msg))


def command(cmd, quiet=False, **kwargs):
    """Echo and execute specified cmd.  If string, execute through BASH.  In that case, cmd could be a pipeline.  Raise an exception if exit code non-zero.  Set check=False if you don't want that exception.  Set quiet=True if you don't want to echo the command.  The result is https://docs.python.org/3/library/subprocess.html#subprocess.CompletedProcess."""
    # This requires version >= 3.6
    assert "shell" not in kwargs, "Please do not override shell.  It is automatically True if command is a string, False if list."
    shell = isinstance(cmd, str)
    if not quiet:
        command_str = cmd if shell else " ".join(cmd)
        tsprint(repr(command_str))
    subproc_args = {
        "check": True
    }
    if shell:
        # By default docker tries a different shell with restricted functionality,
        # and weird handling of pipe errors.  Let's specify BASH here to be safe.
        subproc_args["executable"] = "/bin/bash"
    subproc_args.update(**kwargs)
    return subprocess.run(cmd, shell=shell, **subproc_args)


def command_output(cmd, quiet=False, **kwargs):
    """Echo and execute specified cmd; capture and return its stdout."""
    return command(cmd, quiet=quiet, stdout=subprocess.PIPE, **kwargs).stdout.decode('utf-8')


def backtick(cmd, **kwargs):
    """Execute specified cmd without echoing it.  Capture and return its stdout, stripping any whitespace."""
    return command_output(cmd, quiet=True, **kwargs).strip()


if __name__ == "__main__":
    tsprint(f"Hello from {backtick('pwd')}.  Put tests here.")
