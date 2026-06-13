"""Module providing regression test execution and completion checking."""

import os
import subprocess
import shutil
import re
import shlex
import sys
import time


ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


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

    def __init__(self, test, argv):
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

        if self.argv.engine in ["module", "jupyter"]:
            cmd_exe = sys.executable
        else:
            cmd_exe = f"{self.argv.exe} -i"

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
                else:
                    cmd += arg + " "
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
                                        cwd=test.file_dir,
                                        env=process_env,
                                        shell=True,
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
            self.PerformChecks()
            test.ran = True
            running = False

        return running

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
            with open(stdout_filename, "r", encoding="utf-8") as self.stdout_file:
                shutil.copyfileobj(self.stdout_file, unified_file)
            with open(stderr_filename, "r", encoding="utf-8") as self.stderr_file:
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
            for check in test.checks:
                verbose = self.argv.verbose
                check_passed = check.PerformCheck(output_filename, error_code, verbose)
                passed = passed and check_passed
                if not hasattr(test, "check_results"):
                    test.check_results = []
                test.check_results.append(
                    {"passed": check_passed, "summary": check.GetLastResultSummary()})

                check_annotations = check.GetAnnotations()
                for ann in check_annotations:
                    test.annotations.append(ann)
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
