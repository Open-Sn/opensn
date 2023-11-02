import json
import os
import warnings
from . import checks
import shutil
import time

from . import test_slot


# ========================================================= class definition
class TestConfiguration:
    """Data structure to hold all the necessary info to define a test and
       it's checks"""

    def __init__(self, file_dir: str,
                 filename: str,
                 outfileprefix: str,
                 num_procs: int,
                 checks_params: list,
                 message_prefix: str,
                 dependency: str,
                 args: list,
                 weight_class: str,
                 skip: str):
        """Constructor. Load checks into the data structure"""
        self.file_dir = file_dir
        self.filename = filename
        self.outfileprefix = outfileprefix
        self.num_procs = num_procs
        self.weight_class = weight_class  # default "short"
        self.checks = []
        self.ran = False
        self.submitted = False
        self.annotations = []
        self.dependency = dependency
        self.args = args
        self.skip = skip

        check_num = 0
        for check_params in checks_params:
            check_num += 1
            if not isinstance(check_params, dict):
                warnings.warn(message_prefix + f'Check number {check_num} ' +
                              'is not a dictionary')
                continue

            if "type" not in check_params:
                warnings.warn(message_prefix + f'Check number {check_num} ' +
                              'requires "type" field')
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
        """Shorthand utility get a relative path to a test"""
        return os.path.relpath(self.file_dir + self.filename)

    def GetOutFilenamePrefix(self) -> str:
        if self.outfileprefix == "":
            return self.filename
        else:
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
        """Loops through a test coonfiguration and checks whether a
           dependency has executed"""
        if self.dependency == "":
            return True
        for test in tests:
            if test.filename == self.dependency:
                if test.ran:
                    return True

        return False


# ========================================================= Parse JSON configs
def ParseTestConfiguration(file_path: str):
    """Parses a JSON configuration at the path specified"""
    test_objects = []
    file = open(file_path)
    data = json.load(file)
    file.close()

    err_read = "Error reading " + file_path + ": "

    if not isinstance(data, list):
        warnings.warn(err_read + "Main block is not a list")
        return []

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

        args = []
        if "args" in test_block and not isinstance(test_block["args"], list):
            warnings.warn(message_prefix + '"args" field must be a list')
            continue
        if "args" in test_block:
            args = test_block["args"]

        dependency = ""
        if "dependency" in test_block:
            dependency = test_block["dependency"]

        weight_class = "short"
        if "weight_class" in test_block:
            if isinstance(test_block["weight_class"], str):
                input_weight_class = test_block["weight_class"]
                allowable_list = ["short", "intermediate", "long"]
                if input_weight_class not in allowable_list:
                    warnings.warn(message_prefix + '"weight_class" field, with ' +
                                  f'value "{input_weight_class}" must be in the ' +
                                  'list: ' + allowable_list.__str__())
                    continue
                weight_class = input_weight_class
            else:
                warnings.warn(message_prefix + '"weight_class" field must be a str')
                continue

        outfileprefix = ""
        if "outfileprefix" in test_block:
            if isinstance(test_block["outfileprefix"], str):
                outfileprefix = test_block["outfileprefix"]
            else:
                warnings.warn(message_prefix + '"outfileprefix" field must be a str')
                continue

        skip_reason = ""
        if "skip" in test_block:
            if isinstance(test_block["skip"], str):
                input_reason = test_block["skip"]
                if len(input_reason) == 0:
                    warnings.warn(message_prefix + '"skip" field must be a '
                                                   'non-zero length str')
                    continue
                skip_reason = test_block["skip"]
            else:
                warnings.warn(message_prefix + '"skip" field must be a str')
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
                                         skip=skip_reason)
            test_objects.append(new_test)
        except ValueError:
            continue

    return test_objects


# ========================================================= Top level
#                                                           functions
def ListFilesInDir(folder: str, ext=None):
    """Lists the files in a directory, non-recursively. Optionally an extension
       can be used as a filter"""
    files = []
    dirs_and_files = os.listdir(folder)
    for item in dirs_and_files:
        if not os.path.isdir(item):
            if ext is None:
                files.append(item)
            else:
                base_name, extension = os.path.splitext(item)
                if extension == ext:
                    files.append(item)
    return files


def BuildSearchHierarchyForTests(argv):
    """Finds lua files recursively and creates a map of directories to tests"""
    test_dir = argv.directory

    if not os.path.isdir(test_dir):
        raise Exception('"' + test_dir + '" directory does not exist')

    test_hierarchy = {}  # Map of directories to lua files
    for dir_path, sub_dirs, files in os.walk(test_dir):
        for file_name in files:
            base_name, extension = os.path.splitext(file_name)
            if extension == ".lua":
                abs_dir_path = os.path.abspath(dir_path) + "/"
                if abs_dir_path not in test_hierarchy:
                    test_hierarchy[abs_dir_path] = [file_name]
                else:
                    test_hierarchy[abs_dir_path].append(file_name)

    return test_hierarchy


def ConfigureTests(test_hierarchy: dict, argv):
    """Search through a map of dirs-to-lua-file and looks for a .json file
       that will then be used to create a test object. Also preps the
       out and gold dirs"""

    specific_test = ""
    if argv.test is not None:
        specific_test = argv.test
        print("specific_test=" + specific_test)

    test_objects = []
    for testdir in test_hierarchy:
        for config_file in ListFilesInDir(testdir, ".json"):
            sub_test_objs = ParseTestConfiguration(testdir + config_file)
            for obj in sub_test_objs:
                if specific_test != "" and obj.filename != specific_test:
                    print("skipping " + obj.filename)
                    continue
                test_objects.append(obj)

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
        print(f"test {test_num} " + test.__str__())


def RunTests(tests: list, argv):
    """Actually runs the tests. This routine dynamically checks the system
       load in order to use the system maximally"""
    start_time = time.perf_counter()

    capacity = max(4, argv.jobs)
    system_load = 0

    test_slots = []

    specific_test = ""
    if argv.test is not None:
        specific_test = argv.test

    weight_class_map = ["long", "intermediate", "short"]
    weight_classes_allowed = []
    if 0 <= argv.weights <= 7:
        binary_weights = '{0:03b}'.format(argv.weights)
        for k in range(0, 3):
            if binary_weights[k] == '1':
                weight_classes_allowed.append(weight_class_map[k])
    else:
        warnings.warn('Illegal value "' + str(argv.weights) + '" supplied ' +
                      'for argument -w, --weights')

    print("Executing tests with class in: " + weight_classes_allowed.__str__())

    k = 0
    while True:
        k += 1

        done = True
        # Check for tests to run
        for test in tests:
            if test.ran or (test.weight_class not in weight_classes_allowed):
                continue
            done = False

            if not test.submitted and test.CheckDependencies(tests):
                if test.num_procs <= (capacity - system_load):
                    system_load += test.num_procs

                    new_slot = test_slot.TestSlot(test, argv)

                    # This will only run if a specific test has been
                    # specified
                    if new_slot.test.filename == specific_test:
                        print("Running " + new_slot.test.GetTestPath() + ":")

                    test_slots.append(new_slot)

        # Check test progress
        system_load = 0
        for slot in test_slots:
            if slot.Probe():
                system_load += slot.test.num_procs

        time.sleep(0.01)
        # print(system_load)

        if done:
            break

    num_tests_failed = 0
    for slot in test_slots:
        if not slot.passed:
            num_tests_failed += 1

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time

    print("Done executing tests with class in: " +
          weight_classes_allowed.__str__())

    num_skipped_tests = 0
    for test in tests:
        if not test.ran:
            num_skipped_tests += 1

    if num_skipped_tests > 0:
        print()
        print(f"\033[93mNumber of skipped tests : {num_skipped_tests}\033[0m")

        if num_skipped_tests > 0:
            print("\033[93mSkipped tests:")
            for test in tests:
                if not test.ran:
                    print(test.filename + f' class="{test.weight_class}"')
            print("\033[0m", end='')

    print()
    print("Elapsed time            : {:.2f} seconds".format(elapsed_time))
    print(f"Number of tests run     : {len(test_slots)}")
    print(f"Number of failed tests  : {num_tests_failed}")

    if num_tests_failed > 0:
        return 1
    return 0


def PrintCaughtWarnings(warning_manager, name: str):
    if len(warning_manager) > 0:
        print(f"{name}:")
    for w in warning_manager:
        print("\033[93m" + str(w.category) + "\033[0m", end='')
        print(" ", w.message)
