/**\page DevManTestSystem The Test System

\tableofcontents

Tests are extremely important to the overall vision of the project. They can
be tests to see if a simulation behaves in a certain way, tests to check input
language is handled appropriately, tests to check that operations produce the
expected result, and even checks for certain error or warning behavior. In a
nutshell... it assures that you know if something breaks.

The test system is contained in the `test` directory of the project.
\verbatim
├── doc
├── external
├── framework
├── modules
├── resources
├── test    <-- Here
├── tutorials
├── CMakeLists.txt
├── LICENSE
├── README.md
└── configure.sh
\endverbatim

Within this directory we have the `run_tests` script.
\verbatim
test
├── bin             <-- Test executable in here
├── framework
├── modules
├── src             <-- main.cc in here
├── CMakeLists.txt
└── run_tests       <-- Primary script
\endverbatim

\section DevManTestSystem_sec1 1 Compiling the tests
The test sources, contained throughout the `test` directory, is compiled along
with the regular project but it does not form part of the library (it is a
separate executable). We can do this because ChiTech's compile time is super
short, and ddditionally, the benefit is that we can get compiler errors if we
break any interfaces.

The entry point for the test executable is the `main.cc` contained in the
`test/src` directory. Thereafter all other tests sources are added to the
executable and linked together using \ref DevManStaticRegistration.

The executable is called `opensn-test` and is contained in the `test` directory where OpenSn was
built.

\section DevManTestSystem_sec2 2 What is a test?
Comprises up to 4 things:
- A directory in which the test resides. Multiple tests can be in the same
  directory.
- A `.lua` input file that will initiate the tests.
- A `.json` configuration file specifying one or more tests with associated
  checks. Our convention is to have just one of these in a folder and to
  name it `YTests.json` (the Y always places it at the bottom of the directory)
- Optional. A `.cc` file implementing specific unit testing code.
- Optional. A `.gold` file if the test involves a gold-file check.

Example test:
```
test/example/
├── example_test.lua
├── example_test.cc
└── YTests.json
```
Here we have `example_test.lua` that contains only a single line:
\include test/example/example_test.lua
which executes a wrapped function defined in `example_test.cc`
\include test/example/example_test.cc

The `YTests.json` has the following syntax:
\include test/example/YTests.json

Where we introduce our JSON test configuration option.
- The JSON file starts with a square-bracket wrap `[...]` defining an array of
 <B>Test Blocks</B>.
- Each <B>Test Block</B> is wrapped with curly-braces `{...}`. For this example
 the test had the following parameters:
 - `"comment"` = Free to use to annotate a test. Does not get used internally.
 - `"file"` = The name of the `.lua` that initiates the test.
 - `"num_procs"` = The number of mpi processes to use.
 - `"checks"` = An array of checks `[..Checks..]`, where we have the following
   checks:
   - `StrCompare` check. Checks for the presence of a key string.
   - `ErrorCode` check. Checks for a specified exit code. In this case 0,
     meaning successful execution.

\subsection DevManTestSystem_sec2_1 2.1 Test-blocks documentation
Specifies a specific test
Parameters:
- `"file"` : The name of the `.lua` that initiates the test.
- `"num_procs"` : The number of mpi processes to use.
- `"checks"` : An array of checks `[..Checks..]`
- `"args"` : An array of arguments to pass to the executable
- `"weight_class"` : An optional string, either "short", "intermediate" or
   "long" used to filter different
   tests length. The default is "short". "long" should be used for tests >2min.
- `"outfileprefix"` : Optional parameter. Will default to `"file"` but can be
  to change the output file name (outfileprefix+".out") so that the same input
  file can be used for different tests.
- `"skip"` : Optional parameter. Must be non-empty string stating the reason the
  test was skipped. The presence of this string causes the test to be skipped.

All other keys are ignored so feel free to peruse something like `"comment"` to
annotate the test.

\subsection DevManTestSystem_sec2_2 2.2 Checks documentation
Currently we have the following tests:
- \ref DevManTestSystem_sec2_2_1
- \ref DevManTestSystem_sec2_2_2
- \ref DevManTestSystem_sec2_2_3
- \ref DevManTestSystem_sec2_2_4
- \ref DevManTestSystem_sec2_2_5
- \ref DevManTestSystem_sec2_2_6
\subsubsection DevManTestSystem_sec2_2_1 2.2.1 KeyValuePairCheck
Looks for a key with a floating point value right after it.\n
Parameters:
- `"type"` : "KeyValuePair"
- `"key"` : The key-string to look for
- `"goldvalue"` : Float value
- `"tol"` : Tolerance on the goldvalue
- `"skip_lines_until"` : Optional. Do not check lines in the output file until
  this string is encountered.
  e.g. `"skip_lines_until": "LinearBoltzmann::KEigenvalueSolver execution"`. This
  is useful if a simulation is expected to have multiples of the key-string but
  you only want the last one.

\subsubsection DevManTestSystem_sec2_2_2 2.2.2 StrCompareCheck
Can do one of two things, 1) looks for the presence of the key, returns success
if it is present, or 2) looks for the presence of the key and if present (and the
`"wordnum"` parameter is present then checks if the `"wordnum"`-th word equals
that specified by the `"gold"` value.\n
Parameters:
- `"type"` : "StrCompare"
- `"key"` : The key-string to look for
- `"wordnum"` : Optional. If supplied then "gold" needs to specified.
- `"gold"` : Golden word
- `"skip_lines_until"` : Optional. Do not check lines in the output file until
  this string is encountered.
  e.g. `"skip_lines_until": "LinearBoltzmann::KEigenvalueSolver execution"`. This
  is useful if a simulation is expected to have multiples of the key-string but
  you only want the last one.

\subsubsection DevManTestSystem_sec2_2_3 2.2.3 FloatCompareCheck
On the line containing the key, compares the `"wordnum"`-th word against the
specified gold-value.\n
Parameters:
- `"type"` : "FloatCompare"
- `"key"` : The key-string to look for
- `"wordnum"` : The word number on the line containing the key that will be
                used in the check.
- `"gold"` : Golden value (`float`)
- `"tol"` : The floating point tolerance to use
- `"skip_lines_until"` : Optional. Do not check lines in the output file until
  this string is encountered.
  e.g. `"skip_lines_until": "LinearBoltzmann::KEigenvalueSolver execution"`. This
  is useful if a simulation is expected to have multiples of the key-string but
  you only want the last one.

\subsubsection DevManTestSystem_sec2_2_4 2.2.4 IntCompareCheck
Integer version of \ref DevManTestSystem_sec2_2_3. On the line containing the
key, compares the `"wordnum"`-th word against the
specified gold-value.\n
Parameters:
- `"type"` : "IntCompare"
- `"key"` : The key-string to look for
- `"wordnum"` : The word number on the line containing the key that will be
                used in the check.
- `"gold"` : Golden value (`int`)
- `"skip_lines_until"` : Optional. Do not check lines in the output file until
  this string is encountered.
  e.g. `"skip_lines_until": "LinearBoltzmann::KEigenvalueSolver execution"`. This
  is useful if a simulation is expected to have multiples of the key-string but
  you only want the last one.

\subsubsection DevManTestSystem_sec2_2_5 2.2.5 ErrorCodeCheck
Compares the return/error code of the test with a specified value.\n
Parameters:
- `"type"` : "ErrorCode"
- `"error_code"` : The return code required to pass

\subsubsection DevManTestSystem_sec2_2_6 2.2.6 GoldFileCheck
Compares the contents of the test output to a golden output file.\n
Parameters:
- `"type"` : "GoldFile"
- `"scope_keyword"` : Optional. Restrict the gold comparison to within a section
  of the respective gold/output file that are between the keywords
  `<scope_keyword>_BEGIN` and `<scope_keyword>_END`.
- `"candidate_filename"` : Optional. If supplied, this check will use this file
  rather than the test's output file. For example if `zorba.csv` is provided then
  `zorba.csv` will be compared against `zorba.csv.gold`.
- `"skiplines_top"`: Number of lines at the top of both the gold and comparison
  file to skip in the comparison check.

\section DevManTestSystem_sec3 3 Running the test system
The tests are executed by executing the `run_tests` script, for example:
```
> test/run_tests -j 8 --directory test
```
or
```
> test/run_tests -j 8 -d test
```
Example portion of the output:
```
Description of what you are looking at:
[XX]<test directory>/<lua filename>.lua.......[annotations]Passed/Failed
[XX] = number of mpi processes
[annotations] = error messages preventing a test from running
[lua file missing] = A test pointing to a .lua file was indicated in a .json file but the actual lua file is missing.
[Gold file missing] = Test with a GoldFile check has no .gold file in the gold/ directory. If the input is Input.lua
                      then there needs to be a Input.lua.gold file in the gold directory. If this was the first run
                      then copy Input.lua.out from the out/ directory and use that as the gold.
[Python error] = A python error occurred. Run with -v 1 to see the error.

[ 1]test/framework/chi_misc_utils/chi_misc_utils_test_00.lua......................................................Passed
[ 1]test/framework/parameters/params_test_00.lua..................................................................Passed
[ 1]test/framework/parameters/params_test_01a.lua.................................................................Passed
[ 1]test/framework/parameters/params_test_01b.lua.................................................................Passed
[ 1]test/framework/parameters/params_test_01c.lua.................................................................Passed
[ 2]test/framework/chi_data_types/chi_data_types_test_00.lua......................................................Passed
[ 1]test/framework/parameters/params_test_01d.lua.................................................................Passed
[ 4]test/framework/chi_mesh/ReadWavefrontObj1.lua.................................................................Passed
[ 1]test/framework/tutorials/fv_test1.lua.........................................................................Passed
[ 1]test/framework/tutorials/fv_test2.lua.........................................................................Passed
[ 1]test/framework/tutorials/tutorial_06_wdd.lua..................................................................Passed
[ 4]test/framework/tutorials/pwlc_test1.lua.......................................................................Passed
[ 1]test/framework/chi_math/chi_math_test_00.lua..................................................................Passed

```
The `run_tests` script can be executed on any folder within the `test` directory.
This is a really great feature, it means that you can restrict your testing
to specific areas of the code. For example, if you know you only made changes
to a specific module then there is no need to rerun the framework tests.

The `run_tests` has a number of useful arguments:
\verbatim
options:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        The test directory to process
  -t TEST, --test TEST  A specific test to run
  --exe EXE             The executable to use for testing
  -j JOBS, --jobs JOBS  Allow N jobs at once
  -v VERBOSE, --verbose VERBOSE
                        Controls verbose failure
\endverbatim

The functionality here allows one to execute only a subset of tests. For
example, to only execute the framework tests we can do
```
> test/run_tests -j 8 -d test/framework
```

If interested in a specific test you can narrow the tests even further:
```
> test/run_tests -j 8 -d test/framework -t params_test_00.lua
```
 */
