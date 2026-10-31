"""Module providing regression test execution and completion checking."""

import os
import re
import shlex
import shutil
import signal
import stat
import subprocess
import sys
import time
import uuid

ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")

# Each test runs in a private working directory beside its input's directory, named
# .<directory name>.run-<id>, so the files tests write never collide. The working directory holds
# symbolic links to the input directory's entries, so its other entries are the files the test
# wrote.
RUN_DIRECTORY_RE = re.compile(r"^\..+\.run-[0-9a-f]+$")
# Directories of runner and Python output, which the test tree check skips
OUTPUT_DIRECTORIES = {"out", "__pycache__"}
# Input-directory entries that are not linked into working directories
UNLINKED_ENTRIES = OUTPUT_DIRECTORIES | {"gold", "tests.json"}


def IsRunnerDirectory(name: str) -> bool:
    """Returns True for the names of output and working directories, which the runner owns."""
    return name in OUTPUT_DIRECTORIES or RUN_DIRECTORY_RE.match(name) is not None


def SavedFilesDirectory(test) -> str:
    """Returns out/<output prefix>.files/ in the input's directory, which receives the files a test
    wrote when it failed or when other tests depend on it. The files of a passed test are read-only,
    so a dependent cannot change them through its links. The next run of the test clears it."""
    return os.path.join(test.file_dir, "out", f"{test.GetOutFilenamePrefix()}.files")


def UntrackedEntries(directory: str) -> tuple:
    """Returns the absolute paths of the untracked entries under directory that git ignores, and of
    those it does not ignore, as two sets. A directory with no tracked files is listed instead of
    its contents. Both sets are empty outside a git work tree."""
    directory = os.path.abspath(directory)

    def ListUntracked(*options) -> set:
        try:
            result = subprocess.run(["git", "-C", directory, "ls-files", "-z", "--others",
                                     "--exclude-standard", "--directory", *options],
                                    capture_output=True, text=True, check=False)
        except FileNotFoundError:
            return set()
        if result.returncode != 0:
            return set()
        return {os.path.normpath(os.path.join(directory, path))
                for path in result.stdout.split("\0") if path}

    return ListUntracked("--ignored"), ListUntracked()


def InputEntries(test, ignored: set) -> list:
    """Returns the names of the input-directory entries that the working directory of a test links
    to. Untracked entries that git ignores, such as files left by runs before working directories,
    are not linked, except the test's own input, which the jupyter engine generates."""
    source = test.file_dir.rstrip("/")
    return sorted(entry for entry in os.listdir(source)
                  if entry not in UNLINKED_ENTRIES
                  and (entry == test.filename or os.path.join(source, entry) not in ignored))


def SnapshotTree(directory: str, excluded=frozenset()) -> dict:
    """Returns the size and modification time of every file and directory under directory, except
    the runner's output and working directories, Python caches, and the excluded directories."""
    snapshot = {}
    for dir_path, sub_dirs, files in os.walk(directory):
        sub_dirs[:] = [name for name in sub_dirs if not IsRunnerDirectory(name)
                       and os.path.join(dir_path, name) not in excluded]
        for name in sub_dirs + files:
            path = os.path.join(dir_path, name)
            info = os.lstat(path)
            snapshot[os.path.relpath(path, directory)] = (info.st_size, info.st_mtime_ns)
    return snapshot


def TreeChanges(before: dict, after: dict, directory: str) -> list:
    """Returns the created, removed, and modified paths between two snapshots of directory."""
    return sorted(
        [f"created  {path}" for path in after.keys() - before.keys()]
        + [f"removed  {path}" for path in before.keys() - after.keys()]
        + [f"modified {path}" for path in before.keys() & after.keys()
           if before[path] != after[path] and not os.path.isdir(os.path.join(directory, path))]
    )


def WrittenEntries(run_directory: str) -> list:
    """Returns the paths of the entries a test wrote in its working directory: all entries except
    links and Python caches."""
    paths = [os.path.join(run_directory, entry) for entry in os.listdir(run_directory)
             if entry != "__pycache__"]
    return [path for path in paths if not os.path.islink(path)]


def CreateRunDirectory(test, producer) -> str:
    """Creates the private working directory of a test and returns its path.

    The directory has the same parent as the input's directory, so relative paths to shared assets
    resolve as they do there, and holds symbolic links to the input directory's files (see
    InputEntries). Files a test writes by relative path therefore stay in its own directory. The
    files that the test's dependency (producer) wrote are linked in as well, in place of input
    entries of the same name.
    """
    source = test.file_dir.rstrip("/")
    parent, name = os.path.split(source)
    run_directory = os.path.join(parent, f".{name}.run-{uuid.uuid4().hex[:12]}")
    os.mkdir(run_directory)
    for entry in test.input_entries:
        os.symlink(os.path.join(source, entry), os.path.join(run_directory, entry))
    if producer is not None and os.path.isdir(SavedFilesDirectory(producer)):
        saved_directory = SavedFilesDirectory(producer)
        for entry in os.listdir(saved_directory):
            link = os.path.join(run_directory, entry)
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(os.path.join(saved_directory, entry), link)
    return run_directory


def RemoveRunDirectory(test):
    """Removes the private working directory of a test.

    The files that a failed test, or a test that other tests depend on, wrote are first moved to
    SavedFilesDirectory(test).
    """
    if not test.passed or test.has_dependents:
        written = WrittenEntries(test.run_directory)
        if written:
            saved_directory = SavedFilesDirectory(test)
            os.makedirs(saved_directory, exist_ok=True)
            for path in written:
                shutil.move(path, saved_directory)
            if not test.passed:
                test.saved_output_directory = saved_directory
            else:
                MakeFilesReadOnly(saved_directory)
    shutil.rmtree(test.run_directory)
    test.run_directory = None


def MakeFilesReadOnly(directory: str):
    """Removes write permission from the files under directory. Directories stay writable, so the
    files can still be removed."""
    for dir_path, _, files in os.walk(directory):
        for name in files:
            path = os.path.join(dir_path, name)
            if not os.path.islink(path):
                os.chmod(path, stat.S_IMODE(os.lstat(path).st_mode) & ~0o222)


class InputDirectoryWatch:
    """Detects changes that tests make to their input directories through the links in their working
    directories, for example by writing to an input file or into a linked directory.

    An input directory is compared with its last snapshot whenever one of its tests finishes. A
    change is attributed to every test of that directory that ran since the previous comparison: the
    finished test and those still running. Each of them fails. Input directories nested in this one
    are left to their own comparisons.
    """

    def __init__(self, tests: list, slots: list):
        self.slots = slots
        self.directories = {test.file_dir.rstrip("/") for test in tests}
        self.snapshots = {directory: self._Snapshot(directory) for directory in self.directories}

    def _Snapshot(self, directory: str) -> dict:
        return SnapshotTree(directory, self.directories - {directory})

    def Check(self, test):
        """Compares the input directory of a finished test with its last snapshot."""
        directory = test.file_dir.rstrip("/")
        before = self.snapshots[directory]
        after = self._Snapshot(directory)
        self.snapshots[directory] = after
        changes = TreeChanges(before, after, directory)
        if not changes:
            return
        for slot in self.slots:
            if (slot.test.file_dir == test.file_dir and slot.process is not None
                    and not slot.test.ran):
                slot.test.input_changes += changes


def VisibleLength(text: str) -> int:
    """Returns string length without ANSI color escape sequences."""
    return len(ANSI_RE.sub("", text))


def ShortenMiddle(text: str, max_length: int) -> str:
    """Shortens text to max_length by replacing the middle with an ellipsis."""
    if max_length <= 0 or len(text) <= max_length:
        return text
    if max_length <= 3:
        return text[:max_length]

    keep_left = max(1, (max_length - 3) // 3)
    keep_right = max_length - 3 - keep_left
    return text[:keep_left] + "..." + text[-keep_right:]


def FormatDisplayName(test) -> str:
    """Returns the compact test name used in terminal output."""
    test_path = os.path.join(test.file_dir, test.filename)
    parent_dir = os.path.basename(os.path.dirname(test_path))
    display_name = os.path.join(parent_dir, os.path.basename(test.filename))

    args_str = ', '.join(map(str, test.args))
    if len(args_str) != 0:
        display_name += " (" + args_str + ")"

    if getattr(test, "env", {}):
        env_str = ", ".join(f"{key}={value}" for key, value in sorted(test.env.items()))
        display_name += " {" + env_str + "}"

    return display_name


class TestSlot:
    """Data structure to hold information regarding the execution of a test"""

    def __init__(self, test, argv, watch: InputDirectoryWatch):
        self.watch = watch
        self.process = None
        self.stdout_file = None
        self.stderr_file = None
        self.test = test
        self.passed = False
        self.argv = argv
        self.command: str = ""
        self.start_time = None
        self.end_time = None

        self._Run()

    def _Run(self):
        """Protected method to run the test"""
        test = self.test
        test.submitted = True
        stdout_filename = test.file_dir + f"out/{test.GetOutFilenamePrefix()}.tmp"
        stderr_filename = test.file_dir + f"out/{test.GetOutFilenamePrefix()}.err"

        if test.skip != "":
            return

        test.run_directory = CreateRunDirectory(test, test.producer)

        if self.argv.engine in ["module", "jupyter"]:
            cmd_exe = sys.executable
        else:
            prefix = f"{self.argv.exe_prefix} " if self.argv.exe_prefix else ""
            cmd_exe = f"{prefix}{shlex.quote(self.argv.exe)} -i"

        process_env = os.environ.copy()
        process_env.update(test.env)

        base_name, extension = os.path.splitext(test.filename)
        cmd = self.argv.mpi_cmd + " " + str(test.num_procs) + " "
        cmd += cmd_exe + " "
        cmd += test.filename + " "
        cmd += "--suppress-color "
        if extension == ".py":
            cmd += "--py master_export=False "
            for arg in test.args:
                if arg.find("\"") >= 0:
                    qarg = arg.replace('"', '\\"')
                    if sys.platform.startswith('darwin'):
                        cmd += "--py \"" + qarg + "\" "
                    else:
                        cmd += "--py " + qarg + " "
                elif arg.startswith("-"):
                    # A console option, such as "-v 1"
                    cmd += arg + " "
                else:
                    # A Python assignment, such as "apply_wgdsa=True"
                    cmd += "--py " + arg + " "
        env_prefix = ""
        if test.env:
            env_prefix = " ".join(
                f"{key}={shlex.quote(value)}" for key, value in sorted(test.env.items())
            ) + " "
        self.command = env_prefix + cmd

        # Write command to output file
        self.stdout_file = open(stdout_filename, "w", encoding="utf-8")
        self.stderr_file = open(stderr_filename, "w", encoding="utf-8")

        self.start_time = time.perf_counter()
        self.process = subprocess.Popen(cmd,
                                        cwd=test.run_directory,
                                        env=process_env,
                                        shell=True,
                                        start_new_session=True,
                                        stdout=self.stdout_file,
                                        stderr=self.stderr_file,
                                        universal_newlines=True)

    def Probe(self):
        """Probes the test to see if it finished"""
        test = self.test
        running = True

        if test.ran:
            running = False
            return running

        if not test.ran and test.skip != "":
            self.PerformChecks()
            test.ran = True
            running = False
            return running

        if self.process.poll() is not None:
            self.end_time = time.perf_counter()
            self.CreateUnifiedOutput()
            self.watch.Check(test)
            self.PerformChecks()
            test.ran = True
            running = False
            RemoveRunDirectory(test)

        return running

    def _WaitForGroup(self, timeout: float) -> bool:
        """Waits until no process of the test's group remains. Returns False on timeout."""
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            # Reap the shell so that it does not count as a member
            self.process.poll()
            try:
                os.killpg(self.process.pid, 0)
            except ProcessLookupError:
                return True
            except PermissionError:
                # macOS reports a group whose members are all zombies this way
                return True
            time.sleep(0.05)
        return False

    def Stop(self):
        """Stops the test's process group, including MPI ranks, if it is still running, and waits
        for every process of the group to exit, so that its working directory can be removed"""
        if self.process is None or self.process.poll() is not None:
            return
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(self.process.pid, sig)
            except (ProcessLookupError, PermissionError):
                break
            if self._WaitForGroup(10):
                break
        self.process.wait()
        for file in (self.stdout_file, self.stderr_file):
            if file is not None:
                file.close()

    def CreateUnifiedOutput(self):
        """Combines output and error files into one file"""
        test = self.test
        stdout_filename = test.file_dir + f"out/{test.GetOutFilenamePrefix()}.tmp"
        stderr_filename = test.file_dir + f"out/{test.GetOutFilenamePrefix()}.err"
        unified_filename = test.file_dir + f"out/{test.GetOutFilenamePrefix()}.out"

        # Close stdout and stderr files. This will automatically flush them.
        if self.stdout_file is not None:
            if not self.stdout_file.closed:
                self.stdout_file.close()
        if self.stderr_file is not None:
            if not self.stderr_file.closed:
                self.stderr_file.close()

        # Append stderr and stdout to unified file
        with open(unified_filename, "a", encoding="utf-8") as unified_file:
            unified_file.write(self.command + "\n")
            with open(stdout_filename,
                      "r",
                      encoding="utf-8",
                      errors="backslashreplace") as self.stdout_file:
                shutil.copyfileobj(self.stdout_file, unified_file)
            with open(stderr_filename,
                      "r",
                      encoding="utf-8",
                      errors="backslashreplace") as self.stderr_file:
                shutil.copyfileobj(self.stderr_file, unified_file)
        os.remove(stdout_filename)
        os.remove(stderr_filename)

    def PerformChecks(self):
        """Applies to check-suite for the test"""
        test = self.test
        passed = True
        output_filename = f"{test.file_dir}out/{test.GetOutFilenamePrefix()}.out"
        test.check_results = []

        if test.skip == "":
            error_code = self.process.returncode
            if self.argv.fail_on_nonzero_exit and error_code != 0:
                passed = False
                test.annotations.append(f"Exit code {error_code}")
                test.check_results.append(
                    {
                        "passed": False,
                        "summary": f"Process exit code expected=0 actual={error_code}",
                    }
                )
            for check in test.checks:
                verbose = self.argv.verbose
                check.run_directory = test.run_directory
                check_passed = check.PerformCheck(output_filename, error_code, verbose)
                passed = passed and check_passed
                if not hasattr(test, "check_results"):
                    test.check_results = []
                test.check_results.append(
                    {"passed": check_passed, "summary": check.GetLastResultSummary()})

                check_annotations = check.GetAnnotations()
                for ann in check_annotations:
                    test.annotations.append(ann)
            if test.input_changes:
                passed = False
                test.annotations.append("Changed input directory")
                for change in test.input_changes:
                    test.check_results.append(
                        {"passed": False,
                         "summary": "Input directory changed while the test ran: " + change})
        else:
            test.annotations.append("Skipped")

        test_path = os.path.join(test.file_dir, test.filename)
        test_file_name = FormatDisplayName(test)
        test.display_name = test_file_name

        if not os.path.isfile(test_path):
            test.annotations.append("Input file missing")

        if passed:
            self.passed = True
            status = "\033[32mPassed\033[0m"
        else:
            self.passed = False
            status = "\033[31mFailed\033[0m"

        if test.skip != "":
            status = "\033[36mSkipped\033[0m"

        test.passed = self.passed

        prefix = "\033[33m[{:2d}]\033[0m".format(test.num_procs)

        wall_time_sec = 0.0
        if self.start_time is not None:
            end_time = self.end_time if self.end_time is not None else time.perf_counter()
            wall_time_sec = end_time - self.start_time

        opensn_elapsed_time_sec = None
        if os.path.exists(output_filename):
            with open(output_filename, 'r', encoding='utf-8') as file:
                for line in file:
                    found = re.search("Elapsed execution time:", line)
                    if found:
                        values_slice = re.split(r'[,:]', line.strip())
                        opensn_elapsed_time_sec = (float(values_slice[1]) * 3600
                                                   + float(values_slice[2]) * 60
                                                   + float(values_slice[3]))
                        break

        elapsed_time_sec = opensn_elapsed_time_sec
        if elapsed_time_sec is None:
            elapsed_time_sec = wall_time_sec
        test.elapsed_time_sec = elapsed_time_sec
        time_message = "{:.1f}s".format(elapsed_time_sec)

        annotations = " ".join(f"\033[36m[{annotation}]\033[0m"
                               for annotation in test.annotations)

        terminal_width = min(shutil.get_terminal_size(fallback=(120, 24)).columns, 120)
        status_column = max(40, terminal_width - 26)
        max_name_length = max(10, status_column - VisibleLength(prefix) - 2)
        test_file_name = ShortenMiddle(test_file_name, max_name_length)
        left = prefix + " " + test_file_name
        separator_width = max(1, status_column - VisibleLength(left))
        separator = " " * separator_width
        message = left + separator + status + " " + time_message
        if annotations:
            message += " " + annotations

        print(message)

        if test.skip != "":
            print("\tSkip reason: " + test.skip)
