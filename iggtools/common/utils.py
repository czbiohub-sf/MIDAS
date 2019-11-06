#!/usr/bin/env python3

import os
import sys
import time
import subprocess
import multiprocessing
import random
import traceback
import io
from fnmatch import fnmatch
from functools import wraps


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


class InputStream:
    '''
    Basic usage:

        with InputStream("/path/to/file.lz4") as stream:
            for line in stream:
                tsprint(line)

    To inspect just the first line,

        with InputStream("/path/to/file.lz4") as stream:
            tsprint(stream.readline())
            stream.ignore_errors()

    Failing to read through the entire stream would normally raise an exception,
    unless you call stream.ignore_errors(), as in the above example.

    To run through some custom filters,

        with InputStream("/path/to/file.lz4", """awk '{print $2}' | sed 's=foo=bar='""") as stream:
            for filtered_line in stream:
                tsprint(filtered_line)

    If your filters include commands like "head" or "tail" that can fail to consume all input,
    make sure to call stream.ignore_errors() before you exit the context.

    If the path is in S3, data is streamed from S3 directly, without downloading to
    a local file.  TODO:  In some cases, local caching may be desired.  You can simulate
    that using "tee" in the filter.

    Formats lz4, bz2, gz and plain text are supported.

    Wildcards in path are also supported, but must expand to precisely 1 matching file.
    '''

    def __init__(self, path, filters=None, check_path=True, binary=False):
        if check_path:
            path = smart_glob(path, expected=1)[0]
        cat = 'set -o pipefail; '
        if path.startswith("s3://"):
            cat += f"aws s3 --quiet cp {path} -"
        else:
            cat += f"cat {path}"
        if path.endswith(".lz4"):
            cat += " | lz4 -dc"
        elif path.endswith(".bz2"):
            cat += " | lbzip2 -dc"
        elif path.endswith(".gz"):
            cat += " | gzip -dc"
        if filters:
            cat += f" | {filters}"
        self.cat = cat
        self.path = path
        self.subproc = 0
        self.ignore_called_process_errors = False
        self.binary = binary

    def ignore_errors(self):
        """Exiting the context before consuming all data would normally raise an exception.  To suppress that, call this function."""
        self.ignore_called_process_errors = True

    def __enter__(self):
        self.subproc = command(self.cat, popen=True, stdout=subprocess.PIPE)
        self.subproc.__enter__()
        self.subproc.stdout.ignore_errors = self.ignore_errors
        return self.subproc.stdout if self.binary else text_mode(self.subproc.stdout)  # Note the subject of the WITH statement is the subprocess STDOUT

    def __exit__(self, etype, evalue, etraceback):
        result = self.subproc.__exit__(etype, evalue, etraceback)  # pylint: disable=assignment-from-no-return
        if not self.ignore_called_process_errors:
            returncode = self.subproc.returncode
            if returncode != 0:
                msg = f"Non-zero exit code {returncode} from reader of {self.path}."
                if evalue:
                    # Only a warning as we don't want to silence the other exception.
                    tsprint(f"WARNING: {msg}")
                else:
                    # Surface this as an error.
                    assert returncode == 0, msg
        return result


class OutputStream:
    '''
    Same idea as InputStream, but for output.  Handles compression etc transparently.
    '''

    def __init__(self, path, filters=None):
        if path.startswith("s3://"):
            cat = f"aws s3 --quiet cp - {path}"
        else:
            cat = f"cat > {path}"
        if path.endswith(".lz4"):
            cat = "lz4 -c | " + cat
        elif path.endswith(".bz2"):
            cat = "lbzip2 -c | " + cat
        elif path.endswith(".gz"):
            cat = "gzip -c | " + cat
        if filters:
            cat = f"{filters} | " + cat
        self.cat = 'set -o pipefail; ' + cat
        self.path = path
        self.subproc = 0
        self.ignore_called_process_errors = False

    def ignore_errors(self):
        """Exiting the context before consuming all data would normally raise an exception.  To suppress that, call this function.  This can happen if you specify head or tail as a filter while debugging.  You probably don't want to see this called in production code with OutputStream."""
        self.ignore_called_process_errors = True

    def __enter__(self):
        self.subproc = command(self.cat, popen=True, stdin=subprocess.PIPE)
        self.subproc.__enter__()
        self.subproc.stdin.ignore_errors = self.ignore_errors
        return self.subproc.stdin  # Note the subject of the WITH statement is the subprocess STDIN

    def __exit__(self, etype, evalue, etraceback):
        result = self.subproc.__exit__(etype, evalue, etraceback)  # pylint: disable=assignment-from-no-return
        if not self.ignore_called_process_errors:
            returncode = self.subproc.returncode
            if returncode != 0:
                msg = f"Non-zero exit code {returncode} from reader of {self.path}."
                if evalue:
                    # Only a warning as we don't want to silence the other exception.
                    tsprint(f"WARNING: {msg}")
                else:
                    # Surface this as an error.
                    assert returncode == 0, msg
        return result


def text_mode(stream):
    # Thank you https://stackoverflow.com/questions/31188333/text-mode-adapter-for-binary-or-text-mode-file
    if isinstance(stream, io.TextIOBase):
        result = stream
    elif isinstance(stream, (io.BufferedIOBase, io.RawIOBase)):
        result = io.TextIOWrapper(stream)
    else:
        tsprint(f"WARNING:  Incomplete implementation of utils.text_mode.  Might cause problems if object type {type(stream)} is not already a text-mode stream.")
        result = stream
    return result


def strip_eol(lines_iterable):
    for line in lines_iterable:
        yield line.rstrip('\n')


def command(cmd, quiet=False, popen=False, **kwargs):
    """Echo and execute specified cmd.  If string, execute through BASH.  In that case, cmd could be a pipeline.  Raise an exception if exit code non-zero.  Set check=False if you don't want that exception.  Set quiet=True if you don't want to echo the command.  The result is https://docs.python.org/3/library/subprocess.html#subprocess.CompletedProcess."""
    # This requires version >= 3.6
    assert "shell" not in kwargs, "Please do not override shell.  It is automatically True if command is a string, False if list."
    shell = isinstance(cmd, str)
    if not quiet:
        command_str = cmd if shell else " ".join(cmd)
        tsprint(repr(command_str))
    subproc_args = {}
    if not popen:
        subproc_args["check"] = True
    if shell:
        # By default docker tries a different shell with restricted functionality,
        # and weird handling of pipe errors.  Let's specify BASH here to be safe.
        subproc_args["executable"] = "/bin/bash"
    subproc_args.update(**kwargs)
    if popen:
        return subprocess.Popen(cmd, shell=shell, **subproc_args)
    return subprocess.run(cmd, shell=shell, **subproc_args)


def command_output(cmd, quiet=False, **kwargs):
    """Echo and execute specified cmd; capture and return its stdout."""
    return command(cmd, quiet=quiet, stdout=subprocess.PIPE, **kwargs).stdout.decode('utf-8')


def backtick(cmd, **kwargs):
    """Execute specified cmd without echoing it.  Capture and return its stdout, stripping any whitespace."""
    return command_output(cmd, quiet=True, **kwargs).strip()


def smart_glob(pattern, expected=range(0, sys.maxsize), memory=None):
    """Return list of files that match specified pattern.  Raise exception if item count unexpected.  Memoize if memory dict provided."""
    pdir, file_pattern = pattern.rsplit("/", 1)
    def match_pattern(filename):
        return fnmatch(filename, file_pattern)
    matching_files = list(filter(match_pattern, smart_ls(pdir, memory=memory)))
    actual = len(matching_files)
    expected_str = str(expected)
    if isinstance(expected, int):
        expected = [expected]
    assert actual in expected, f"Expected {expected_str} file(s) for {pattern}, found {actual}"
    return [f"{pdir}/{mf}" for mf in sorted(matching_files)]


find_files = smart_glob


def smart_ls(pdir, missing_ok=True, memory=None):
    "Return a list of files in pdir.  This pdir can be local or in s3.  If memory dict provided, use it to memoize.  If missing_ok=True, swallow errors (default)."
    result = memory.get(pdir) if memory else None
    if result == None:
        try:
            if pdir.startswith("s3://"):
                s3_dir = pdir
                if not s3_dir.endswith("/"):
                    s3_dir += "/"
                output = backtick(["aws", "s3", "ls", s3_dir])
                rows = output.strip().split('\n')
                result = [r.split()[-1] for r in rows]
            else:
                output = backtick(["ls", pdir])
                result = output.strip().split('\n')
        except Exception as e:
            msg = f"Could not read directory: {pdir}"
            if missing_ok and isinstance(e, subprocess.CalledProcessError):
                tsprint(f"INFO: {msg}")
                result = []
            else:
                tsprint(f"ERROR: {msg}")
                raise
        if memory != None:
            memory[pdir] = result
    return result


def parse_table(lines, columns):
    lines = strip_eol(lines)
    headers = next(lines).split('\t')  # pylint: disable=stop-iteration-return
    column_indexes = []
    for c in columns:
        assert c in headers, f"Column not found in table headers: {c}"
        column_indexes.append(headers.index(c))
    for i, l in enumerate(lines):
        values = l.split('\t')
        if len(headers) != len(values):
            assert False, f"Line {i} has {len(values)} columns;  was expecting {len(headers)}."
        yield [values[i] for i in column_indexes]


def retry(operation, MAX_TRIES=3):
    # Note the use of a separate random generator for retries so transient
    # errors won't perturb other random streams used in the application.
    invocation = [0]
    @wraps(operation)
    def wrapped_operation(*args, **kwargs):
        randgen = None
        remaining_attempts = MAX_TRIES
        delay = 2.0
        while remaining_attempts > 1:
            try:
                return operation(*args, **kwargs)
            except:
                if randgen == None:
                    invocation[0] += 1
                    randgen = random.Random(os.getpid() * 10000 + invocation[0]).random
                wait_time = delay * (1.0 + 2.0 * randgen())
                tsprint(f"Sleeping {wait_time} seconds before retry {MAX_TRIES - remaining_attempts + 1} of {operation} with {args}, {kwargs}.")
                time.sleep(wait_time)
                delay *= 3.0
                remaining_attempts -= 1
        # The last attempt is outside try/catch so caller can handle exception
        return operation(*args, **kwargs)
    return wrapped_operation


def suppress_exceptions(func):
    @wraps(func)
    def func_noexc(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            tsprint(f"\nSuppressing exception below and returning None from {func} called with arguments {args} {kwargs}:\n{traceback.format_exc(limit=-1)}")
            return None
    return func_noexc


def portion(iterable, portion_size):
    assert portion_size > 0
    p = []
    for item in iterable:
        p.append(item)
        if len(p) == portion_size:
            yield p
            p = []
    if p:
        yield p


def multi_map(func, items, num_procs):
    p = multiprocessing.Pool(num_procs)
    items = list(items)
    results = p.map(func, items)
    return dict(zip(items, results))


if __name__ == "__main__":
    tsprint(f"Hello from {backtick('pwd')}.  Put tests here.")
