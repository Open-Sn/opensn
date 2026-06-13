"""Module providing job management for regression tests."""

import json
import os
import fnmatch
import warnings
import shutil
import time
from . import checks
from . import test_slot


class TestConfiguration:
    """Data structure that holds the necessary info to define a test and its checks"""

    def __init__(self, file_dir: str,
                 filename: str,
                 outfileprefix: str,
                 num_procs: int,
                 checks_params: list,
                 message_prefix: str,
                 dependency: str,
                 args: list,
                 weight_class: str,
                 skip: str,
                 requires_gpu: bool,
                 env: dict):
        """Constructor. Load checks into the data structure"""
        self.file_dir = file_dir
        self.filename = filename
        self.outfileprefix = outfileprefix
        self.num_procs = num_procs
        self.weight_class = weight_class  # default "short"
        self.checks = []
        self.ran = False
        self.submitted = False
        self.passed = False
        self.annotations = []
        self.dependency = dependency
        self.args = args
        self.skip = skip
        self.requires_gpu = requires_gpu
        self.env = env
        self.skip_reason_override = None

        check_num = 0
        for check_params in checks_params:
            check_num += 1
            if not isinstance(check_params, dict):
                warnings.warn(message_prefix + f'Check number {check_num} ' + 'is not a dictionary')
                continue

            if "type" not in check_params:
                warnings.warn(
                    message_prefix + f'Check number {check_num} ' + 'requires "type" field')
                continue

            if check_params["type"] == "KeyValuePair":
                try:
                    prefix = message_prefix + f'Check number {check_num} '
                    new_check = checks.KeyValuePairCheck(check_params, prefix)
                    self.checks.append(new_check)
                except ValueError:
                    continue
            elif check_params["type"] == "StrCompare":
                try:
                    prefix = message_prefix + f'Check number {check_num} '
                    new_check = checks.StrCompareCheck(check_params, prefix)
                    self.checks.append(new_check)
                except ValueError:
                    continue
            elif check_params["type"] == "FloatCompare":
                try:
                    prefix = message_prefix + f'Check number {check_num} '
                    new_check = checks.FloatCompareCheck(check_params, prefix)
                    self.checks.append(new_check)
                except ValueError:
                    continue
            elif check_params["type"] == "IntCompareCheck":
                try:
                    prefix = message_prefix + f'Check number {check_num} '
                    new_check = checks.IntCompareCheck(check_params, prefix)
                    self.checks.append(new_check)
                except ValueError:
                    continue
            elif check_params["type"] == "ErrorCode":
                try:
                    prefix = message_prefix + f'Check number {check_num} '
                    new_check = checks.ErrorCodeCheck(check_params, prefix)
                    self.checks.append(new_check)
                except ValueError:
                    continue
            elif check_params["type"] == "GoldFile":
                try:
                    prefix = message_prefix + f'Check number {check_num} '
                    new_check = checks.GoldFileCheck(check_params, prefix)
                    self.checks.append(new_check)
                except ValueError:
                    continue
            else:
                warnings.warn("Unsupported check type: " + check_params["type"])
                raise ValueError

        if len(self.checks) == 0:
            warnings.warn(message_prefix + " has no valid checks")
            raise ValueError

    def GetTestPath(self):
        """Shorthand utility to get the relative path to a test"""
        return os.path.relpath(self.file_dir + self.filename)

    def GetOutFilenamePrefix(self) -> str:
        """Get the output filename prefix"""
        if self.outfileprefix == "":
            return self.filename
        return self.outfileprefix

    def __str__(self):
        """Converts the class to a readable format"""
        output = f'file_dir="{self.file_dir}" '
        output += f'filename="{self.filename}" '
        output += f'num_procs={self.num_procs} '

        check_num = 0
        for check in self.checks:
            check_num += 1
            output += f'\n    check_num={check_num} '
            output += check.__str__()

        return output

    def CheckDependencies(self, tests):
        """Loops through a test configuration and checks whether a dependency has executed"""
        if self.dependency is None:
            return True
        for test in tests:
            if test.filename == self.dependency:
                if test.ran and test.passed:
                    return True

        return False


def ResolveDependencyMatches(test_objects: dict, dependency: str) -> list:
    """Returns configured tests whose filename matches the dependency string."""
    return [test for test in test_objects.values() if test.filename == dependency]


# Parse JSON configs
def ParseTestConfiguration(file_path: str):
    """Parses a JSON configuration at the path specified"""
    test_objects = {}

    with open(file_path, 'r', encoding='utf-8') as file:
        data = json.load(file)

    err_read = "Error reading " + file_path + ": "

    if not isinstance(data, list):
        warnings.warn(err_read + "Main block is not a list")
        return {}

    test_num = 0
    for test_block in data:
        test_num += 1
        message_prefix = err_read + f'Test {test_num} '

        if not isinstance(test_block, dict):
            warnings.warn(message_prefix + 'is not a dictionary')
            continue
        if "file" not in test_block:
            warnings.warn(message_prefix + 'does not have key "file"')
            continue
        if "num_procs" not in test_block:
            warnings.warn(message_prefix + 'does not have key "num_procs"')
            continue
        if "checks" not in test_block:
            warnings.warn(message_prefix + 'does not have key "checks"')
            continue

        if not isinstance(test_block["file"], str):
            warnings.warn(message_prefix + '"file" field must be a string')
            continue

        if not isinstance(test_block["checks"], list):
            warnings.warn(message_prefix + '"checks" field must be a list')
            continue

        if not isinstance(test_block["num_procs"], int):
            warnings.warn(message_prefix + '"num_procs" field must be an integer')
            continue
        if test_block["num_procs"] <= 0:
            warnings.warn(message_prefix + '"num_procs" field must be positive')
            continue

        args = []
        if "args" in test_block and not isinstance(test_block["args"], list):
            warnings.warn(message_prefix + '"args" field must be a list')
            continue
        if "args" in test_block:
            args = test_block["args"]

        env = {}
        if "env" in test_block:
            if not isinstance(test_block["env"], dict):
                warnings.warn(message_prefix + '"env" field must be a dictionary')
                continue

            normalized_env = {}
            for key, value in test_block["env"].items():
                if not isinstance(key, str):
                    warnings.warn(message_prefix + '"env" keys must be strings')
                    normalized_env = None
                    break
                if not isinstance(value, (str, int, float, bool)):
                    warnings.warn(
                        message_prefix + f'"env" value for "{key}" must be string-, int-, '
                        + 'float-, or bool-like')
                    normalized_env = None
                    break
                normalized_env[key] = str(value)

            if normalized_env is None:
                continue
            env = normalized_env

        dependency = None
        if "dependency" in test_block:
            if isinstance(test_block["dependency"], str):
                dependency = test_block["dependency"]
            else:
                warnings.warn(message_prefix + '"dependency" field must be a string')
                continue

        weight_class = "short"
        if "weight_class" in test_block:
            if isinstance(test_block["weight_class"], str):
                input_weight_class = test_block["weight_class"]
                allowable_list = ["short", "intermediate", "long"]
                if input_weight_class not in allowable_list:
                    warnings.warn(message_prefix + '"weight_class" field, with '
                                  + f'value "{input_weight_class}" must be in the '
                                  + 'list: ' + str(allowable_list))
                    continue
                weight_class = input_weight_class
            else:
                warnings.warn(message_prefix + '"weight_class" field must be a string')
                continue

        outfileprefix = ""
        if "outfileprefix" in test_block:
            if isinstance(test_block["outfileprefix"], str):
                outfileprefix = test_block["outfileprefix"]
            else:
                warnings.warn(message_prefix + '"outfileprefix" field must be a string')
                continue

        skip_reason = ""
        if "skip" in test_block:
            if isinstance(test_block["skip"], str):
                input_reason = test_block["skip"]
                if len(input_reason) == 0:
                    warnings.warn(message_prefix + '"skip" field must be a non-zero length string')
                    continue
                skip_reason = test_block["skip"]
            else:
                warnings.warn(message_prefix + '"skip" field must be a string')
                continue

        requires_gpu = False
        if "requires_gpu" in test_block:
            gpu_value = test_block["requires_gpu"]
            if isinstance(gpu_value, bool):
                requires_gpu = gpu_value
            elif isinstance(gpu_value, str):
                normalized = gpu_value.strip().lower()
                if normalized in ("true", "1", "yes", "on"):
                    requires_gpu = True
                elif normalized in ("false", "0", "no", "off", ""):
                    requires_gpu = False
                else:
                    warnings.warn(message_prefix + '"requires_gpu" field must be boolean-like')
                    continue
            elif isinstance(gpu_value, (int, float)):
                requires_gpu = bool(gpu_value)
            else:
                warnings.warn(message_prefix + '"requires_gpu" field must be boolean-like')
                continue

        try:
            new_test = TestConfiguration(file_dir=os.path.dirname(file_path) + "/",
                                         filename=test_block["file"],
                                         outfileprefix=outfileprefix,
                                         num_procs=test_block["num_procs"],
                                         checks_params=test_block["checks"],
                                         message_prefix=message_prefix,
                                         dependency=dependency,
                                         args=args,
                                         weight_class=weight_class,
                                         skip=skip_reason,
                                         requires_gpu=requires_gpu,
                                         env=env)
            test_key = (
                new_test.filename,
                tuple(map(str, new_test.args)),
                tuple(sorted(new_test.env.items())),
            )
            if test_key in test_objects:
                warnings.warn(message_prefix + "duplicates an existing file/args test entry")
                continue
            test_objects[test_key] = new_test
        except ValueError:
            continue

    invalid_test_keys = []
    for test_key, test in test_objects.items():
        if test.dependency is None:
            continue

        dependency_matches = ResolveDependencyMatches(test_objects, test.dependency)
        if len(dependency_matches) == 0:
            warnings.warn(err_read + f'test "{test.filename}" depends on missing test '
                          + f'"{test.dependency}"')
            invalid_test_keys.append(test_key)
            continue
        if len(dependency_matches) > 1:
            warnings.warn(err_read + f'test "{test.filename}" has ambiguous dependency '
                          + f'"{test.dependency}". Dependencies must identify exactly one test '
                          + 'entry within the same tests.json file.')
            invalid_test_keys.append(test_key)

    for test_key in invalid_test_keys:
        test_objects.pop(test_key, None)

    return test_objects


def ListFilesInDir(folder: str, ext=None):
    """Lists the files in a directory, non-recursively. A file extension can be used as a filter"""
    files = []
    dirs_and_files = os.listdir(folder)
    for item in dirs_and_files:
        item_path = os.path.join(folder, item)
        if not os.path.isdir(item_path):
            if ext is None:
                files.append(item)
            else:
                base_name, extension = os.path.splitext(item)
                if extension == ext:
                    files.append(item)
    return sorted(files)


def ConvertNbToScript(notebook_path: str) -> str:
    """
    Converts a Jupyter notebook (.ipynb) to a Python script (.py) in the same directory.

    Parameters
    ----------
    notebook_path : str
        Path to the input .ipynb file.

    Returns
    -------
    str
        Path to the generated .py file.
    """
    # Return old file name if it is not a notebook
    if not notebook_path.endswith(".ipynb"):
        return notebook_path

    try:
        from nbconvert import PythonExporter
        import nbformat
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Jupyter notebook tests require nbconvert and nbformat. "
            "Install them or run the regression suite with --engine console."
        ) from exc

    # Read notebook
    nb_node = nbformat.read(notebook_path, as_version=4)
    # Create exporter
    exporter = PythonExporter()
    source_code, _ = exporter.from_notebook_node(nb_node)
    # Determine output .py path
    base, _ = os.path.splitext(notebook_path)
    output_py_path = base + ".py"
    # Save the Python script
    with open(output_py_path, "w") as f:
        f.write(source_code)
    return output_py_path


def BuildSearchHierarchyForTests(argv):
    """Finds input files recursively and creates a map of directories to tests"""
    test_dir = argv.directory

    if not os.path.isdir(test_dir):
        raise Exception('"' + test_dir + '" directory does not exist')

    test_hierarchy = {}  # Map of directories to input files
    for dir_path, sub_dirs, files in os.walk(test_dir):
        sub_dirs.sort()
        files.sort()
        for file_name in files:
            if argv.engine == "jupyter":
                file_name = ConvertNbToScript(os.path.join(dir_path, file_name))
            base_name, extension = os.path.splitext(file_name)
            if extension == ".py":
                abs_dir_path = os.path.abspath(dir_path) + "/"
                if abs_dir_path not in test_hierarchy:
                    test_hierarchy[abs_dir_path] = [file_name]
                else:
                    test_hierarchy[abs_dir_path].append(file_name)

    return test_hierarchy


def SelectTests(test_obj, selector: str, base_directory: str) -> bool:
    """Returns true if a configured test matches a -t/--test selector."""
    test_path = os.path.join(test_obj.file_dir, test_obj.filename)
    rel_path = os.path.relpath(test_path, base_directory)
    candidates = [test_obj.filename, rel_path, test_path]
    return any(
        candidate == selector or fnmatch.fnmatch(candidate, selector) for candidate in candidates
    )


def ConfigureTests(test_hierarchy: dict, argv):
    """Search through a map of dirs-to-input-file and looks for a .json file that will then be used
       to create a test object. Also preps the out and gold directories."""

    specific_tests = argv.test
    if specific_tests is not None:
        print("Selected tests: " + ", ".join(specific_tests))

    test_objects = []
    for testdir in test_hierarchy:
        for config_file in ListFilesInDir(testdir, ".json"):
            sub_test_objs = ParseTestConfiguration(testdir + config_file)
            specific_test_dependencies = []
            for obj in sub_test_objs.values():
                if argv.gpu:
                    if not obj.requires_gpu:
                        continue
                else:
                    if obj.requires_gpu:
                        continue
                if specific_tests is None or any(
                    SelectTests(obj, selector, argv.directory) for selector in specific_tests
                ):
                    test_objects.append(obj)
                    if specific_tests is not None and obj.dependency is not None:
                        specific_test_dependencies.append(obj.dependency)
                elif argv.show_filtered:
                    print("skipping " + obj.filename)

            # If a specific test has dependencies, also add them to the list of executed tests
            for dependency in specific_test_dependencies:
                dependency_matches = ResolveDependencyMatches(sub_test_objs, dependency)
                if len(dependency_matches) == 0:
                    warnings.warn("Specified dependency '" + dependency + "' does not exist.")
                    continue
                if len(dependency_matches) > 1:
                    warnings.warn("Specified dependency '" + dependency + "' is ambiguous.")
                    continue

                dependency_obj = dependency_matches[0]
                include_dependency = (
                    dependency_obj.requires_gpu if argv.gpu else not dependency_obj.requires_gpu
                )
                if include_dependency and dependency_obj not in test_objects:
                    test_objects.append(dependency_obj)
                elif not include_dependency:
                    warnings.warn("Dependency '" + dependency + "' filtered by GPU mode.")

        # If the out directory exists then we clear it
        if os.path.isdir(testdir + "out/"):
            shutil.rmtree(testdir + "out/")

        # If the out directory does not exist then we create it
        if not os.path.isdir(testdir + "out/"):
            os.mkdir(testdir + "out/")

        # If the gold directory does not exist then we create it
        if not os.path.isdir(testdir + "gold/"):
            os.mkdir(testdir + "gold/")

    return test_objects


def EchoTests(tests: list):
    """For debugging, echos the string format of each test configuration"""
    test_num = 0
    for test in tests:
        print(f"test {test_num} " + str(test))


def RunTests(tests: list, argv):
    """Actually runs the tests. This routine dynamically checks the system load."""
    start_time = time.perf_counter()

    # os.cpu_count() may not be ideal in this case since it returns the number of logical cpus.
    capacity = os.cpu_count()
    if argv.jobs > 0:
        capacity = argv.jobs
    capacity = max(1, capacity)
    system_load = 0
    test_slots = []
    scheduler_blocked = False

    weight_class_map = ["long", "intermediate", "short"]
    weight_classes_allowed = []
    if 0 <= argv.weights <= 7:
        binary_weights = '{0:03b}'.format(argv.weights)
        for k in range(0, 3):
            if binary_weights[k] == '1':
                weight_classes_allowed.append(weight_class_map[k])
    else:
        warnings.warn('Illegal value "' + str(argv.weights) + '" supplied '
                      + 'for argument -w, --weights')

    runnable_tests = [
        test for test in tests
        if test.weight_class in weight_classes_allowed
        and ((test.requires_gpu and argv.gpu) or (not test.requires_gpu and not argv.gpu))
    ]
    max_test_procs = max((test.num_procs for test in runnable_tests), default=1)
    if max_test_procs > capacity:
        print(
            "\033[93mNote:\033[0m requested job capacity "
            f"({capacity}) is smaller than the largest selected test "
            f"({max_test_procs} MPI processes). Raising capacity to {max_test_procs}."
        )
        capacity = max_test_procs

    print("Executing tests")
    print("  weights  : " + ", ".join(weight_classes_allowed))
    print(f"  selected : {len(runnable_tests)}")
    print(f"  capacity : {capacity} MPI processes")
    print()

    while True:
        # Check test progress
        system_load = 0
        for slot in test_slots:
            if slot.Probe():
                system_load += slot.test.num_procs

        remaining_tests = [
            test for test in tests
            if not test.ran
            and test.weight_class in weight_classes_allowed
            and ((test.requires_gpu and argv.gpu) or (not test.requires_gpu and not argv.gpu))
        ]

        if len(remaining_tests) == 0:
            break

        launched_test = False
        blocked_reasons = {}
        for test in remaining_tests:
            if test.submitted:
                blocked_reasons[test] = "scheduler bookkeeping error after test submission"
                continue

            if not test.CheckDependencies(tests):
                blocked_reasons[test] = f'dependency "{test.dependency}" did not pass'
                continue

            if test.num_procs > (capacity - system_load):
                blocked_reasons[test] = (
                    f"waiting for {test.num_procs} MPI processes; "
                    f"{capacity - system_load} currently available"
                )
                continue

            system_load += test.num_procs
            new_slot = test_slot.TestSlot(test, argv)
            test_slots.append(new_slot)
            launched_test = True

        if system_load == 0 and not launched_test:
            scheduler_blocked = True
            warnings.warn(
                "Stopping test scheduler because remaining tests cannot be launched. "
                "This indicates either an unsatisfied dependency or an internal scheduler error."
            )
            for test in remaining_tests:
                reason = blocked_reasons.get(
                    test, "scheduler could not determine the blocking reason"
                )
                test.skip_reason_override = reason
                warnings.warn("Blocked test: " + GetDisplayName(test))
                warnings.warn("  Reason: " + reason)
            break
        time.sleep(0.1)

    num_tests_failed = 0
    for slot in test_slots:
        if not slot.passed:
            num_tests_failed += 1

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time

    print()
    print("Done executing tests")

    PrintSkippedSummary(tests, weight_classes_allowed, argv)

    print()
    print("Summary")
    print("  elapsed time       : {:.2f} seconds".format(elapsed_time))
    print(f"  tests run          : {len(test_slots)}")
    print(f"  tests failed       : {num_tests_failed}")

    PrintFailureSummary(test_slots)
    PrintSlowestTests(test_slots, argv.slowest)

    if num_tests_failed > 0 or scheduler_blocked:
        return 1
    return 0


def FormatElapsedTime(seconds: float) -> str:
    """Formats elapsed seconds for compact terminal summaries."""
    if seconds < 60.0:
        return f"{seconds:.1f}s"
    if seconds < 3600.0:
        return f"{seconds / 60.0:.1f}m"
    return f"{seconds / 3600.0:.1f}h"


def GetDisplayName(test) -> str:
    """Returns a compact display name for summaries."""
    if hasattr(test, "display_name"):
        return test.display_name
    return test_slot.FormatDisplayName(test)


def GetSkippedReason(test, weight_classes_allowed: list, argv) -> str:
    """Returns why a configured test did not run."""
    if test.skip_reason_override is not None:
        return test.skip_reason_override
    if test.weight_class not in weight_classes_allowed:
        return f'weight "{test.weight_class}" not enabled'
    if test.requires_gpu and not argv.gpu:
        return "requires GPU"
    if argv.gpu and not test.requires_gpu:
        return "non-GPU test filtered by --gpu"
    if test.dependency is not None:
        return f'dependency "{test.dependency}" was not run'
    return "not scheduled"


def PrintSkippedSummary(tests: list, weight_classes_allowed: list, argv):
    """Prints skipped tests grouped by reason."""
    skipped_tests = [test for test in tests if not test.ran]
    if len(skipped_tests) == 0:
        return

    grouped_skips = {}
    for test in skipped_tests:
        reason = GetSkippedReason(test, weight_classes_allowed, argv)
        grouped_skips.setdefault(reason, []).append(test)

    print()
    print(f"\033[93mSkipped tests: {len(skipped_tests)}\033[0m")
    for reason, reason_tests in sorted(grouped_skips.items()):
        print(f"  {len(reason_tests):4d}  {reason}")
        if argv.verbose:
            for test in reason_tests:
                print(f"        {GetDisplayName(test)}")


def PrintFailureSummary(test_slots: list):
    """Prints failed tests and failed check summaries."""
    failed_slots = [slot for slot in test_slots if not slot.passed]
    if len(failed_slots) == 0:
        return

    print()
    print("\033[31mFailed tests:\033[0m")
    for slot in failed_slots:
        test = slot.test
        output_filename = os.path.join(test.file_dir, "out", f"{test.GetOutFilenamePrefix()}.out")
        print(f"  {GetDisplayName(test)}")
        print(f"    output: {output_filename}")

        check_results = getattr(test, "check_results", [])
        failed_checks = [result for result in check_results if not result["passed"]]
        for result in failed_checks:
            summary = result["summary"] if result["summary"] else "<no summary>"
            print(f"    check : {summary}")

        if len(failed_checks) == 0 and test.annotations:
            print("    annotations: " + ", ".join(test.annotations))


def PrintSlowestTests(test_slots: list, max_tests: int = 10):
    """Prints the slowest completed tests."""
    if max_tests <= 0:
        return

    timed_slots = [
        slot for slot in test_slots
        if hasattr(slot.test, "elapsed_time_sec") and slot.test.elapsed_time_sec is not None
    ]
    if len(timed_slots) == 0:
        return

    timed_slots.sort(key=lambda slot: slot.test.elapsed_time_sec, reverse=True)
    timed_slots = timed_slots[:max_tests]

    print()
    print(f"Slowest tests ({len(timed_slots)}):")
    for slot in timed_slots:
        elapsed = FormatElapsedTime(slot.test.elapsed_time_sec)
        print(f"  {elapsed:>8}  {GetDisplayName(slot.test)}")


def PrintCaughtWarnings(warning_manager, name: str):
    """Prints warnings"""
    if len(warning_manager) > 0:
        print(f"{name}:")
    for w in warning_manager:
        print("\033[93m" + str(w.category) + "\033[0m", end='')
        print(" ", w.message)
