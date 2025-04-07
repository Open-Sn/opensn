"""Module providing regression test execution and completion checking."""

import os
import subprocess
import shutil
import re
import sys


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

        self._Run()

    def _Run(self):
        """Protected method to run the test"""
        test = self.test
        test.submitted = True
        stdout_filename = test.file_dir + f"out/{test.GetOutFilenamePrefix()}.tmp"
        stderr_filename = test.file_dir + f"out/{test.GetOutFilenamePrefix()}.err"

        if test.skip != "":
            return

        base_name, extension = os.path.splitext(test.filename)
        cmd = self.argv.mpi_cmd + " " + str(test.num_procs) + " "
        cmd += self.argv.exe + " "
        cmd += "-i " + test.filename + " "
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
        elif extension == ".lua":
            cmd += "--lua master_export=false "
            for arg in test.args:
                if arg.find("\"") >= 0:
                    qarg = arg.replace('"', '\\"')
                    if sys.platform.startswith('darwin'):
                        cmd += "--lua " + "\\\"" + qarg + "\\\" "
                    else:
                        cmd += "--lua " + qarg + " "
                else:
                    cmd += arg + " "
        self.command = cmd

        # Write command to output file
        self.stdout_file = open(stdout_filename, "w", encoding="utf-8")
        self.stderr_file = open(stderr_filename, "w", encoding="utf-8")

        self.process = subprocess.Popen(cmd,
                                        cwd=test.file_dir,
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

        if test.skip == "":
            error_code = self.process.returncode
            for check in test.checks:
                verbose = self.argv.verbose
                check_passed = check.PerformCheck(output_filename, error_code, verbose)
                passed = passed and check_passed

                check_annotations = check.GetAnnotations()
                for ann in check_annotations:
                    test.annotations.append(ann)
        else:
            test.annotations.append("Skipped")

        test_path = os.path.join(test.file_dir, test.filename)
        test_file_name = os.path.relpath(test_path, self.argv.directory)

        if not os.path.isfile(test_path):
            test.annotations.append("Lua file missing")

        args_str = ', '.join(map(str, test.args))
        if len(args_str) != 0:
            args_str = " (" + args_str + ")"
        test_file_name += args_str

        pad = 0
        if passed:
            self.passed = True
            message = "\033[32mPassed\033[0m"
            pad += 5 + 4
        else:
            self.passed = False
            message = "\033[31mFailed\033[0m"
            pad += 5 + 4

        prefix = "\033[33m[{:2d}]\033[0m".format(test.num_procs)
        pad += 5 + 4

        for annotation in test.annotations:
            message = f"\033[36m[{annotation}]\033[0m" + " " + message
            pad += 5 + 4

        width = 120 - len(prefix + test_file_name) + pad
        message = message.rjust(width, ".")

        opensn_elapsed_time_sec = 0.0
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

        time_message = " {:.1f}s".format(opensn_elapsed_time_sec)

        print(prefix + " " + test_file_name + message + time_message)

        if test.skip != "":
            print("\tSkip reason: " + test.skip)
