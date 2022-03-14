#!/usr/bin/env python3

import os
import sys
import time
import subprocess
import json
import multiprocessing
from multiprocessing.pool import ThreadPool
import random
import traceback
import io
from fnmatch import fnmatch
from functools import wraps

# Thread-safe and timestamped prints.
tslock = multiprocessing.RLock()

num_vcpu = multiprocessing.cpu_count()
num_physical_cores = (num_vcpu + 1) // 2


COMPRESSORS = {
    ".lz4": "lz4 -dc",
    ".bz2": "bzip2 -dc",   # consider lbzip2 if available
    ".gz": "gzip -dc"
}


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
        sys.stderr.flush()


def tsprint(msg):
    tserr(tsfmt(msg))


class InputStream:
    '''
    Basic usage:

        with InputStream("/path/to/file.lz4") as stream:
            for line in stream:
                tsprint(line)

    To run through some custom filters,

        with InputStream("/path/to/file.lz4", """awk '{print $2}' | sed 's=foo=bar='""") as stream:
            for filtered_line in stream:
                tsprint(filtered_line)

    If your filters include commands like "head" or "tail" that intentionally skip some input,
    make sure to call stream.ignore_errors() before you exit the context, because, without
    that, failing to read through the entire stream would raise an exception.

    For instance, to inspect just the first line,

        with InputStream("/path/to/file.lz4") as stream:
            tsprint(stream.readline())
            stream.ignore_errors()

    or, similarly, to inspect just line 125,

        with InputStream("/path/to/file.lz4", "head -125 | tail -1") as stream:
            tsprint(stream.readline())
            stream.ignore_errors()

    If the path is in S3, data is streamed from S3 directly, without downloading to
    a local file.  TODO:  In some cases, local caching may be desired.  You can simulate
    that using "tee" in the filter.

    Formats lz4, bz2, gz and plain text are supported.

    Wildcards in path are also supported, but must expand to precisely 1 matching file.
    '''

    def __init__(self, path, through=None, filters=None, check_path=True, binary=False):
        if through != None:
            assert filters == None
            filters = through
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
        self.stream = None

    def ignore_errors(self):
        """Exiting the context before consuming all data would normally raise an exception.  To suppress that, call this function."""
        self.ignore_called_process_errors = True

    def __enter__(self):
        self.subproc = command(self.cat, popen=True, stdout=subprocess.PIPE)
        self.subproc.__enter__()
        self.stream = self.subproc.stdout if self.binary else text_mode(self.subproc.stdout)  # Note the subject of the WITH statement is the subprocess STDOUT
        self.stream.ignore_errors = self.ignore_errors
        return self.stream

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

    def __init__(self, path, through=None, filters=None, binary=False):
        if through != None:
            assert filters == None
            filters = through
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
        self.binary = binary
        self.stream = None

    def ignore_errors(self):
        """Exiting the context before consuming all data would normally raise an exception.  To suppress that, call this function.  This can happen if you specify head or tail as a filter while debugging.  You probably don't want to see this called in production code with OutputStream."""
        self.ignore_called_process_errors = True

    def __enter__(self):
        self.subproc = command(self.cat, popen=True, stdin=subprocess.PIPE)
        self.subproc.__enter__()
        self.stream = self.subproc.stdin if self.binary else text_mode(self.subproc.stdin)  # Note the subject of the WITH statement is the subprocess STDIN
        self.stream.ignore_errors = self.ignore_errors
        return self.stream

    def __exit__(self, etype, evalue, etraceback):
        self.stream.flush()
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


def command(cmd, quiet=True, popen=False, **kwargs):
    """Echo and execute specified cmd.  If string, execute through BASH.  In that case, cmd could be a pipeline.  Raise an exception if exit code non-zero.  Set check=False if you don't want that exception.  Set quiet=True if you don't want to echo the command.  The result is https://docs.python.org/3/library/subprocess.html#subprocess.CompletedProcess."""
    # This requires version >= 3.6
    assert "shell" not in kwargs, "Please do not override shell.  It is automatically True if command is a string, False if list."
    shell = isinstance(cmd, str)
    if not quiet:
        command_str = cmd if shell else " ".join(cmd)
        hack_remove_prefix = "set -o pipefail; "
        print_str = command_str
        if print_str.startswith(hack_remove_prefix):
            print_str = print_str[len(hack_remove_prefix):]
        tsprint(repr(print_str))
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
    quiet = kwargs.pop('quiet', True)
    return command_output(cmd, quiet=quiet, **kwargs).strip()


def smart_glob(pattern, expected=range(0, sys.maxsize), memory=None):
    """Return list of files that match specified pattern.  Raise exception if item count unexpected.  Memoize if memory dict provided."""
    if "/" not in pattern:
        pattern = "./" + pattern
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
                # tsprint(f"INFO: {msg}.  It is OK for that directory to be missing.")
                result = []
            else:
                tsprint(f"ERROR: {msg}")
                raise
        if memory != None:
            memory[pdir] = result
    return result


def select_from_tsv(input_lines, selected_columns=None, schema=None, result_structure=tuple, default_column_type=str):
    """
    Parse selected_columns from a table represented in tab-separated-values format (TSV).

    For example, given an input_stream of lines

        "name\tage\theight\n",
        "tommy\t9\t120\n",
        "masha\t7\t122\n",
        ...

    the generator

        select_from_tsv(input_stream, selected_columns=["name", "height"])

    would yield

       ("tommy", "120"),
       ("masha", "122"),
       ...

    Note how in the above example all values returned are strings.

    Column type information may be specified as follows

        select_from_tsv(input_stream, selected_columns={"name": str, "height": float})

    yielding

        ("tommy", 120.0),
        ("masha", 122.0),
        ...

    If the first line in the input does not list all column headers, the schema must be provided as an additional argument

        select_from_tsv(input_stream, selected_columns=["name", "height"], schema={"name": str, "age": int, "height": float})

    Specifyng a schema argument means the first line of input should contain values, not columnd headers.

    Type information may be specified in either the selected_columns or the schema argument.

    Requesting result_structure=dict will change the output to

       {"name": "tommy", "height": 120.0},
       {"name": "masha", "height": 122.0},
       ...

    When a schema argument is provided, leaving selected_columns unspecified is equivalent to selecting all columns from the schema.

    At least one of the schema or selected_columns arguments must be specified, even if the schema could be inferred from the input headers.  This makes applications more robust.
    """
    assert schema != None or selected_columns != None, "at least one of these arguments must be specified"
    lines = strip_eol(input_lines)
    j = 0
    if schema == None:
        # the first line is expected to list all column names (headers)
        headers = next(lines).split('\t')  # pylint: disable=stop-iteration-return
        schema = headers
        j = 1
    if not isinstance(schema, dict):
        schema = {c: default_column_type for c in schema}
    headers = list(schema.keys())
    # Ensure "selected_columns" is an ordered dict of {column_name: column_type} and "schema" is a superdict of "columns".
    if selected_columns == None:
        # Return all columns
        selected_columns = dict(schema)
    if not isinstance(selected_columns, dict):
        selected_columns = {c: default_column_type for c in selected_columns}
    # Merge column type information from schema and selected_columns
    for c in schema:
        if schema[c] == default_column_type:
            schema[c] = selected_columns.get(c, schema[c])
    for c in selected_columns:
        if selected_columns[c] == default_column_type:
            selected_columns[c] = schema.get(c, selected_columns[c])
    column_indexes = []
    for c in selected_columns:
        assert c in headers, f"Column not found in table headers: {c}"
        column_indexes.append(headers.index(c))
        assert schema[c] == selected_columns[c], f"Conflicting types for column {c}."
    column_types = list(selected_columns.values())
    column_names = list(selected_columns.keys())
    for i, l in enumerate(lines):
        values = l.split('\t')
        if len(headers) != len(values):
            assert False, f"Line {i + j} has {len(values)} columns;  was expecting {len(headers)}."
        # type-convert and reorder values as specified in schema
        ordered_values = (ctype(values[ci]) for ci, ctype in zip(column_indexes, column_types))
        if result_structure in (tuple, list):
            yield result_structure(ordered_values)
        else: # dict
            yield result_structure((c, val) for c, val in zip(column_names, ordered_values))


def _test_select_from_tsv():
    # Unit tests for the select_from_tsv function
    unlabeled_test_data = [
        "tommy\t9\t120\n",
        "masha\t7\t122\n",
    ]
    test_data_with_headers = [
        "name\tage\theight"
    ] + unlabeled_test_data

    result_tuples = [
        ("tommy", 120.0),
        ("masha", 122.0)
    ]
    result_tuples_height_name = [
        (120.0, "tommy"),
        (122.0, "masha")
    ]
    result_dicts_age_name = [
        {"age": 9, "name": "tommy"},
        {"age": 7, "name": "masha"}
    ]

    r = list(select_from_tsv(unlabeled_test_data, ["name", "height"], {"name": str, "age": int, "height": float}))
    assert r == result_tuples, json.dumps(r, indent=4)

    r = list(select_from_tsv(unlabeled_test_data, {"name": str, "height": float}, ["name", "age", "height"]))
    assert r == result_tuples, json.dumps(r, indent=4)

    r = list(select_from_tsv(test_data_with_headers, {"height": float, "name": str}))
    assert r == result_tuples_height_name, json.dumps(r, indent=4)

    r = list(select_from_tsv(test_data_with_headers, {"age": int, "name": str}, result_structure=dict))
    assert r == result_dicts_age_name, json.dumps(r, indent=4)


# Just run this every time the module is imported.   Eventually put in the test section in main below.
_test_select_from_tsv()


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
                t_start = time.time()
                return operation(*args, **kwargs)
            except:
                t_end = time.time()
                if randgen == None:
                    invocation[0] += 1
                    randgen = random.Random(os.getpid() * 10000 + invocation[0]).random
                if t_end - t_start > 30:
                    # For longer operations, the delay should increase, so that the backoff will meaningfully reduce load on the failing service.
                    delay = (t_end - t_start) * 0.2
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


def split(iterable, portion_size):
    assert portion_size > 0
    p = []
    for item in iterable:
        p.append(item)
        if len(p) == portion_size:
            yield p
            p = []
    if p:
        yield p


# The next few functions present a uniform interface to a parallel version of map
# using either multiple threads (preferred) or multiple processes (necessary when
# the function's performance is bottlenecked on CPU-intensive python code).
#
# Use multiprocessing *only* if the Python code itself is slow and bottlenecked on CPU.
# If the problem is just I/O latency, use multithreading, rather than multiprocessing.
# If the CPU intensive work can be done in subcommands, use multithreading rather than
# multiprocessing, to invoke multiple subcommands in parallel.
#
# When multiprocessing, it pays little to exceed the number of physical cores,
# and never more than 2x, particularly if there are merging costs.
#
# With multithreading, on the other hand, there is no reason to constrain
# the number of threads;  by default we set it equal to the number of items.


def hashmap(func, items):
    items = list(items)
    results = map(func, items)
    return dict(zip(items, results))


# private! use multiprocessing_map or multithreading_map instead
def _multi_map(func, items, num_procs, PoolClass):
    p = PoolClass(num_procs)
    return p.map(func, items, chunksize=1)


# private! use multiprocessing_hashmap or multithreading_hashmap instead
def _multi_hashmap(func, items, num_procs, PoolClass):
    items = list(items)
    results = _multi_map(func, items, num_procs, PoolClass)
    return dict(zip(items, results))


# use this *only* if func is CPU bound
def multiprocessing_map(func, items, num_procs=num_physical_cores):
    return _multi_map(func, items, num_procs, multiprocessing.Pool)


# use this if func is not CPU bound
def multithreading_map(func, items, num_threads=None):
    if not num_threads:
        items = list(items)
        num_threads = max(1, len(items))
    return _multi_map(func, items, num_threads, ThreadPool)


# use this *only* if func is CPU bound
def multiprocessing_hashmap(func, items, num_procs=num_physical_cores):
    return _multi_hashmap(func, items, num_procs, multiprocessing.Pool)


# use this if func is not CPU bound
def multithreading_hashmap(func, items, num_threads=None):
    if not num_threads:
        items = list(items)
        num_threads = max(1, len(items))
    return _multi_hashmap(func, items, num_threads, ThreadPool)


def transpose(list_of_tuples):
    # zip is its own inverse, for small enough data
    # this converts [(a, 1), (b, 2), (c, 3)] into ([a, b, c], [1, 2, 3])
    return zip(*list_of_tuples)


def sorted_dict(d):
    return {k: d[k] for k in sorted(d.keys())}


def reordered_dict(d, key_order):
    return {k: d[k] for k in key_order}


def flatten(l):
    return [item for sublist in l for item in sublist]


@retry
def upload(src, dst, check=True):
    assert os.path.isfile(src)
    try:
        command(f"set -o pipefail; lz4 -c {src} | aws s3 cp --only-show-errors - {dst}", check=check)
    except:
        command(f"aws s3 rm {dst}", check=False)
        raise


def upload_star(srcdst):
    src, dst = srcdst
    return upload(src, dst)


def pythonpath():
    # Path from which this program can be called with "python3 -m midas2"
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def datecode(t=None, local=False):
    """The date code incorporates both the unix time and the date in either the local timezone (respecting DST) or UTC (aka GMT).  It is convenient for naming ops subdirectories or tagging events for easy sorting."""
    if t == None:
        t = time.time()
    zt = time.localtime(t) if local else time.gmtime(t)
    return "{year}_{month:02d}_{day:02d}__{t}".format(
        year=zt.tm_year,
        month=zt.tm_mon,
        day=zt.tm_mday,
        t=int(t)
    )


class TimedSection:

    def __init__(self, section_name, quiet=True):
        self.section_name = section_name
        self.t_start = None
        self.quiet = quiet

    def __enter__(self):
        self.t_start = time.time()
        if not self.quiet:
            tsprint(f"Entering section '{self.section_name}'")
        return self

    def __exit__(self, _etype, _evalue, _etraceback):
        tsprint(f"Took {(time.time() - self.t_start)/60:.1f} minutes for {self.section_name}")
        return False  # Do not suppress exceptions


def uncompressed(filename):
    fn, ext = os.path.splitext(filename)
    if ext in COMPRESSORS:
        return fn, COMPRESSORS[ext]
    return filename, "cat"


# Occasional failures in aws s3 cp require a retry.
@retry
def download_reference(ref_path, local_dir="."):
    local_path = os.path.join(local_dir, os.path.basename(ref_path))
    local_path, uncompress_cmd = uncompressed(local_path)
    if os.path.exists(local_path):
        # TODO:  Reuse instead of re-download.  Requires versioning.
        tsprint(f"Overwriting pre-existing {local_path} with reference download.")
        command(f"rm -f {local_path}")
    if not os.path.exists(local_dir):
        command(f"mkdir -p {local_dir}")
    try:
        command(f"set -o pipefail; aws s3 cp --only-show-errors --no-sign-request {ref_path} - | {uncompress_cmd} > {local_path}")
    except:
        command(f"rm -f {local_path}")
        raise
    return local_path


def cat_files(files_of_chunks, one_file, number_of_chunks=20):
    for temp_files in split(files_of_chunks, number_of_chunks):
        command("cat " + " ".join(temp_files) + f" >> {one_file}", quiet=True)


def drop_lz4(filename):
    assert filename.endswith(".lz4")
    return filename[:-4]

# -------------- testing testing testing ---------------------

def random_sleep(args):
    t = time.time()
    idx, quiet = args
    delay = 0.1 + random.random() * 0.3
    if not quiet:
        tsprint(f"{idx:5} falling asleep for {delay:.3f} seconds.")
    time.sleep(delay)
    elapsed = time.time() - t
    overslept = elapsed - delay
    if not quiet and overslept >= 0.1 * delay:
        tsprint(f"{idx:5} overslept by {overslept:.3f} seconds.")
    return (t, idx, elapsed)

def test_sleepers(n, quiet=False):
    tsprint(f"Testing multiprocesing.map order with {n} sleepers on {num_physical_cores} physical cores.")
    t = time.time()
    info = multiprocessing_map(random_sleep, ((idx, quiet) for idx in range(n)))
    order = [idx for _, idx, _ in sorted(info)]
    elapsed_sum = sum(elapsed for _, _, elapsed in info)
    elapsed = time.time() - t
    if not quiet:
        tsprint(f"Elapsed {elapsed:.1f} seconds with effective parallelism {elapsed_sum / max(elapsed, 1e-9):.1f}x.")
    if order != list(range(n)):
        examples = dict([(i, order[i]) for i in range(n) if i != order[i]][:10])
        assert False, f"Error: Out-of-order dispatch: examples {examples}: full order {order}."
    tsprint(f"SUCCESS")

if __name__ == "__main__":
    tsprint(f"Hello from {backtick('pwd')}.  Put tests here.")
    try:
        n_sleepers = int(sys.argv[1])
    except:
        n_sleepers = num_physical_cores * 31
    for n in range(max(1, n_sleepers - 3), n_sleepers + 4):
        test_sleepers(n, quiet=True)
