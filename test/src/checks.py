"""Module providing regression output checking."""

import os.path
import warnings
import re
import pathlib
import difflib
import numbers


def CheckUnknownFields(params: dict, allowed_fields: set, message_prefix: str):
    """Raises when a check definition contains unsupported fields."""
    unknown_fields = set(params.keys()) - allowed_fields
    if len(unknown_fields) > 0:
        warnings.warn(message_prefix + 'Unsupported field(s): '
                      + ', '.join(sorted(unknown_fields)))
        raise ValueError


def IsRealNumber(value) -> bool:
    """Returns true for JSON integer/float values, but not booleans."""
    return isinstance(value, numbers.Real) and not isinstance(value, bool)


def RequireRealNumber(params: dict, field: str, message_prefix: str):
    """Validates that a required field is a real number."""
    if field not in params:
        warnings.warn(message_prefix + f'Missing "{field}" field')
        raise ValueError
    if not IsRealNumber(params[field]):
        warnings.warn(message_prefix + f'"{field}" field must be numeric')
        raise ValueError


def RequireInteger(params: dict, field: str, message_prefix: str):
    """Validates that a required field is an integer."""
    if field not in params:
        warnings.warn(message_prefix + f'Missing "{field}" field')
        raise ValueError
    if not isinstance(params[field], int) or isinstance(params[field], bool):
        warnings.warn(message_prefix + f'"{field}" field must be an integer')
        raise ValueError


def SplitOutputWords(text: str) -> list:
    """Splits output text on whitespace, commas, and equals signs."""
    return re.split(r'\s+|,+|=+', text.rstrip())


def SplitPostKeyWords(text: str) -> list:
    """Splits text after a key and removes separator-created empty fields."""
    return [word for word in SplitOutputWords(text.strip()) if word != ""]


class Check:
    """Base class for a Check data structure"""

    def __init__(self):
        self.annotations = []
        self.last_result_summary = ""

    def __str__(self):
        return "Check base class"

    def PerformCheck(self, filename, errorcode, verbose: bool):
        warnings.warn("Unimplemented base class routine")
        return False

    def GetAnnotations(self):
        return self.annotations

    def SetLastResultSummary(self, summary: str):
        self.last_result_summary = summary

    def GetLastResultSummary(self):
        return self.last_result_summary


class KeyValuePairCheck(Check):
    """Given a string key, checks the floating point value of the word immediately following the
    key"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.key: str = ""
        self.goldvalue: float = 0.0
        self.rel_tol: float = 0.
        self.abs_tol: float = 0.
        self.skip_lines_until: str = ""

        CheckUnknownFields(params,
                           {"type", "key", "goldvalue", "rel_tol",
                            "abs_tol", "skip_lines_until"},
                           message_prefix)
        if "key" not in params:
            warnings.warn(message_prefix + 'Missing "key" field')
            raise ValueError
        if not isinstance(params["key"], str):
            warnings.warn(message_prefix + '"key" field must be a string')
            raise ValueError
        RequireRealNumber(params, "goldvalue", message_prefix)
        if "rel_tol" not in params and "abs_tol" not in params:
            warnings.warn(message_prefix + 'Must specify "rel_tol" and/or "abs_tol" field')
            raise ValueError

        self.key = params["key"]
        self.goldvalue = params["goldvalue"]
        if "abs_tol" in params:
            if not IsRealNumber(params["abs_tol"]):
                warnings.warn(message_prefix + '"abs_tol" field must be numeric')
                raise ValueError
            self.abs_tol = params["abs_tol"]
        else:
            self.abs_tol = 0.
        if "rel_tol" in params:
            if not IsRealNumber(params["rel_tol"]):
                warnings.warn(message_prefix + '"rel_tol" field must be numeric')
                raise ValueError
            self.rel_tol = params["rel_tol"]
        else:
            self.rel_tol = 0.

        if "skip_lines_until" in params:
            val = params["skip_lines_until"]
            if not isinstance(val, str):
                raise ValueError("Expected 'skip_lines_until' value to be str")
            self.skip_lines_until = val

    def __str__(self):
        return 'type="KeyValuePair", ' + f'key="{self.key}", ' + \
            f'goldvalue={self.goldvalue}, ' + \
            f'rel_tol={self.rel_tol}, ' + f'abs_tol={self.abs_tol}'

    def PerformCheck(self, filename, errorcode, verbose: bool):
        self.SetLastResultSummary(f'KeyValuePair key="{self.key}" result=FAILED reason=unknown')
        try:
            with open(filename, "r", encoding='utf-8') as file:
                skip_lines = False
                perform_skip_check = bool(self.skip_lines_until != "")
                if perform_skip_check:
                    skip_lines = True

                found_key = False
                last_value = None
                last_line = ""
                lines = file.readlines()
                for line in lines:
                    if perform_skip_check:
                        if line.find(self.skip_lines_until) >= 0:
                            skip_lines = False
                            perform_skip_check = False
                    if skip_lines:
                        continue

                    key_pos = line.find(self.key)
                    if key_pos >= 0:
                        found_key = True
                        postkey = line[(key_pos + len(self.key)):].strip()

                        words = SplitPostKeyWords(postkey)
                        if len(words) == 0:
                            raise ValueError("word split failure: " + line)

                        try:
                            value = float(words[0])
                        except Exception:
                            self.annotations.append("Python error")
                            self.SetLastResultSummary(
                                f'KeyValuePair key="{self.key}" expected={self.goldvalue:.9g} '
                                f'result=FAILED reason=parse_error')
                            if verbose:
                                warnings.warn('Failed to convert word "' + words[0] + '" to float\n'
                                              + 'post key: ' + postkey + "\n"
                                              + 'words' + str(words))

                            return False

                        last_value = value
                        last_line = line
                        if abs(value - self.goldvalue) <= self.abs_tol + self.rel_tol * \
                           abs(self.goldvalue):
                            self.SetLastResultSummary(
                                f'KeyValuePair key="{self.key}" expected={self.goldvalue:.9g} '
                                f'actual={value:.9g} result=PASSED')
                            return True
                        self.SetLastResultSummary(
                            f'KeyValuePair key="{self.key}" expected={self.goldvalue:.9g} '
                            f'actual={value:.9g} abs_tol={self.abs_tol:.3g} '
                            f'rel_tol={self.rel_tol:.3g} result=FAILED')

                if verbose and found_key:
                    print("Check failed : " + self.__str__() + "\n" + last_line)
                elif verbose:
                    print('Check failed : key, "' + self.key + '", not found ')
                if not found_key:
                    self.SetLastResultSummary(
                        f'KeyValuePair key="{self.key}" expected={self.goldvalue:.9g} '
                        f'result=FAILED reason=key_not_found')
                elif last_value is not None:
                    self.SetLastResultSummary(
                        f'KeyValuePair key="{self.key}" expected={self.goldvalue:.9g} '
                        f'actual={last_value:.9g} abs_tol={self.abs_tol:.3g} '
                        f'rel_tol={self.rel_tol:.3g} result=FAILED')
        except FileNotFoundError as e:
            warnings.warn(str(e))
            self.SetLastResultSummary(
                f'KeyValuePair key="{self.key}" expected={self.goldvalue:.9g} '
                f'result=FAILED reason=file_not_found')
        except Exception as e:
            self.annotations.append("Python error")
            self.SetLastResultSummary(
                f'KeyValuePair key="{self.key}" expected={self.goldvalue:.9g} '
                f'result=FAILED reason=python_error')
            if verbose:
                warnings.warn(str(e))

        return False


class StrCompareCheck(Check):
    """Given a key to identify a line, compares the specific word number against
       a golden value. The word is expected to be a string"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.key: str = ""
        self.wordnum: int = -1
        self.gold: str = ""
        self.skip_lines_until: str = ""

        CheckUnknownFields(params, {"type", "key", "wordnum", "gold", "skip_lines_until"},
                           message_prefix)
        if "key" not in params:
            warnings.warn(message_prefix + 'Missing "key" field')
            raise ValueError
        if not isinstance(params["key"], str):
            warnings.warn(message_prefix + '"key" field must be a string')
            raise ValueError

        if "wordnum" in params:
            RequireInteger(params, "wordnum", message_prefix)
            if "gold" not in params:
                warnings.warn(message_prefix + 'Missing "gold" field')
                raise ValueError
            if not isinstance(params["gold"], str):
                warnings.warn(message_prefix + '"gold" field must be a string')
                raise ValueError

        self.key = params["key"]
        if "wordnum" in params:
            self.wordnum = params["wordnum"]
            self.gold = params["gold"]

        if "skip_lines_until" in params:
            val = params["skip_lines_until"]
            if not isinstance(val, str):
                raise ValueError("Expected 'skip_lines_until' value to be str")
            self.skip_lines_until = val

    def __str__(self):
        return 'type="StrCompare", ' + f'key="{self.key}", ' + \
            f'wordnum={self.wordnum}, ' + \
            f'gold={self.gold} '

    def PerformCheck(self, filename, errorcode, verbose: bool):
        self.SetLastResultSummary(f'StrCompare key="{self.key}" result=FAILED reason=unknown')
        try:
            with open(filename, "r", encoding='utf-8') as file:
                skip_lines = False
                perform_skip_check = bool(self.skip_lines_until != "")
                if perform_skip_check:
                    skip_lines = True

                lines = file.readlines()
                for line in lines:
                    if perform_skip_check:
                        if line.find(self.skip_lines_until) >= 0:
                            skip_lines = False
                            perform_skip_check = False
                    if skip_lines:
                        continue

                    key_pos = line.find(self.key)
                    if key_pos >= 0:
                        if self.wordnum < 0:
                            self.SetLastResultSummary(
                                f'StrCompare key="{self.key}" result=PASSED found=true')
                            return True
                        words = SplitOutputWords(line)

                        if len(words) <= self.wordnum:
                            warnings.warn("word count: " + str(len(words)) + "\n"
                                          + "line: " + line.rstrip() + "\n"
                                          + "words: " + str(words))
                            raise ValueError(f"Required word {self.wordnum} does not exist")

                        value = words[self.wordnum]

                        if value == self.gold:
                            self.SetLastResultSummary(
                                f'StrCompare key="{self.key}" expected="{self.gold}" '
                                f'actual="{value}" result=PASSED')
                            return True
                        elif verbose:
                            warnings.warn("Check failed : " + self.__str__() + "\n"
                                          + line + "\n"
                                          + str(words))
                        self.SetLastResultSummary(
                            f'StrCompare key="{self.key}" expected="{self.gold}" '
                            f'actual="{value}" result=FAILED')

                if verbose:
                    warnings.warn('Check failed : key, "' + self.key + '", not found '
                                  + str(self.wordnum))
                self.SetLastResultSummary(
                    f'StrCompare key="{self.key}" expected="{self.gold}" '
                    f'result=FAILED reason=key_not_found')
        except FileNotFoundError as e:
            print(str(e))
            self.SetLastResultSummary(
                f'StrCompare key="{self.key}" expected="{self.gold}" '
                f'result=FAILED reason=file_not_found')
        except Exception as e:
            self.annotations.append("Python error")
            self.SetLastResultSummary(
                f'StrCompare key="{self.key}" expected="{self.gold}" '
                f'result=FAILED reason=python_error')
            if verbose:
                print(str(e))

        return False


class FloatCompareCheck(Check):
    """Given a key to identify a line, compares the specific word number against
       a golden value. The word is expected to be a float"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.key: str = ""
        self.wordnum: int = 0
        self.gold: float = 0.0
        self.abs_tol: float = 0.
        self.rel_tol: float = 0.
        self.skip_lines_until: str = ""

        CheckUnknownFields(params,
                           {"type", "key", "wordnum", "gold", "rel_tol",
                            "abs_tol", "skip_lines_until"},
                           message_prefix)
        if "key" not in params:
            warnings.warn(message_prefix + 'Missing "key" field')
            raise ValueError
        if not isinstance(params["key"], str):
            warnings.warn(message_prefix + '"key" field must be a string')
            raise ValueError
        RequireInteger(params, "wordnum", message_prefix)
        RequireRealNumber(params, "gold", message_prefix)
        if "rel_tol" not in params and "abs_tol" not in params:
            warnings.warn(message_prefix + 'Must specify "rel_tol" and/or "abs_tol" field')
            raise ValueError

        self.key = params["key"]
        self.wordnum = params["wordnum"]
        self.gold = params["gold"]
        if "abs_tol" in params:
            if not IsRealNumber(params["abs_tol"]):
                warnings.warn(message_prefix + '"abs_tol" field must be numeric')
                raise ValueError
            self.abs_tol = params["abs_tol"]
        else:
            self.abs_tol = 0.
        if "rel_tol" in params:
            if not IsRealNumber(params["rel_tol"]):
                warnings.warn(message_prefix + '"rel_tol" field must be numeric')
                raise ValueError
            self.rel_tol = params["rel_tol"]
        else:
            self.rel_tol = 0.

        if "skip_lines_until" in params:
            val = params["skip_lines_until"]
            if not isinstance(val, str):
                raise ValueError("Expected 'skip_lines_until' value to be str")
            self.skip_lines_until = val

    def __str__(self):
        return 'type="FloatCompare", ' + f'key="{self.key}", ' + \
            f'wordnum={self.wordnum}, ' + \
            f'gold={self.gold}, ' + \
            f'rel_tol={self.rel_tol}, ' + f'abs_tol={self.abs_tol}'

    def PerformCheck(self, filename, errorcode, verbose: bool):
        self.SetLastResultSummary(f'FloatCompare key="{self.key}" result=FAILED reason=unknown')
        try:
            with open(filename, "r", encoding='utf-8') as file:
                skip_lines = False
                perform_skip_check = bool(self.skip_lines_until != "")
                if perform_skip_check:
                    skip_lines = True

                lines = file.readlines()
                for line in lines:
                    if perform_skip_check:
                        if line.find(self.skip_lines_until) >= 0:
                            skip_lines = False
                            perform_skip_check = False
                    if skip_lines:
                        continue

                    key_pos = line.find(self.key)
                    if key_pos >= 0:
                        words = SplitOutputWords(line)

                        if len(words) <= self.wordnum:
                            warnings.warn("word count: " + str(len(words)) + "\n"
                                          + "line: " + line + "\n"
                                          + "\n"
                                          + "words: " + str(words))
                            raise ValueError(f"Required word {self.wordnum} does not exist")

                        value = 0.0
                        try:
                            value = float(words[self.wordnum])
                        except Exception:
                            self.annotations.append("Python error")
                            self.SetLastResultSummary(
                                f'FloatCompare key="{self.key}" expected={self.gold:.9g} '
                                f'result=FAILED reason=parse_error')
                            if verbose:
                                warnings.warn("    Check failed :\n"
                                              + "    Check = " + self.__str__() + "\n"
                                              + "    line = " + line
                                              + "    words = " + str(words) + "\n"
                                              + "    info = Failed to convert word "
                                              + str(self.wordnum) + " to float.")
                            return False

                        if abs(value - self.gold) <= self.abs_tol + self.rel_tol * abs(self.gold):
                            self.SetLastResultSummary(
                                f'FloatCompare key="{self.key}" expected={self.gold:.9g} '
                                f'actual={value:.9g} result=PASSED')
                            return True
                        elif verbose:
                            warnings.warn("Check failed : " + self.__str__() + "\n"
                                          + line + "\n"
                                          + str(words))
                        self.SetLastResultSummary(
                            f'FloatCompare key="{self.key}" expected={self.gold:.9g} '
                            f'actual={value:.9g} abs_tol={self.abs_tol:.3g} '
                            f'rel_tol={self.rel_tol:.3g} result=FAILED')

                if verbose:
                    warnings.warn('Check failed : key, "' + self.key + '", not found '
                                  + str(self.wordnum))
                self.SetLastResultSummary(
                    f'FloatCompare key="{self.key}" expected={self.gold:.9g} '
                    f'result=FAILED reason=key_not_found')
        except FileNotFoundError as e:
            print(str(e))
            self.SetLastResultSummary(
                f'FloatCompare key="{self.key}" expected={self.gold:.9g} '
                f'result=FAILED reason=file_not_found')
        except Exception as e:
            self.annotations.append("Python error")
            self.SetLastResultSummary(
                f'FloatCompare key="{self.key}" expected={self.gold:.9g} '
                f'result=FAILED reason=python_error')
            if verbose:
                print(str(e))

        return False


class IntCompareCheck(Check):
    """Given a key to identify a line, compares the specific word number against
       a golden value. The word is expected to be a int"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.key: str = ""
        self.wordnum: int = 0
        self.gold: int = 0
        self.skip_lines_until: str = ""

        CheckUnknownFields(params, {"type", "key", "wordnum", "gold", "skip_lines_until"},
                           message_prefix)
        if "key" not in params:
            warnings.warn(message_prefix + 'Missing "key" field')
            raise ValueError
        if not isinstance(params["key"], str):
            warnings.warn(message_prefix + '"key" field must be a string')
            raise ValueError
        RequireInteger(params, "wordnum", message_prefix)
        RequireInteger(params, "gold", message_prefix)

        self.key = params["key"]
        self.wordnum = params["wordnum"]
        self.gold = params["gold"]

        if "skip_lines_until" in params:
            val = params["skip_lines_until"]
            if not isinstance(val, str):
                raise ValueError("Expected 'skip_lines_until' value to be str")
            self.skip_lines_until = val

    def __str__(self):
        return 'type="IntCompare", ' + f'key="{self.key}", ' + \
            f'wordnum={self.wordnum}, ' + \
            f'gold={self.gold} '

    def PerformCheck(self, filename, errorcode, verbose: bool):
        self.SetLastResultSummary(f'IntCompare key="{self.key}" result=FAILED reason=unknown')
        try:
            with open(filename, "r", encoding='utf-8') as file:
                skip_lines = False
                perform_skip_check = bool(self.skip_lines_until != "")
                if perform_skip_check:
                    skip_lines = True

                lines = file.readlines()
                for line in lines:
                    if perform_skip_check:
                        if line.find(self.skip_lines_until) >= 0:
                            skip_lines = False
                            perform_skip_check = False
                    if skip_lines:
                        continue

                    key_pos = line.find(self.key)
                    if key_pos >= 0:
                        words = SplitOutputWords(line)

                        if len(words) <= self.wordnum:
                            warnings.warn("word count: " + str(len(words)) + "\n"
                                          + "line: " + line + "\n"
                                          + "\n"
                                          + "words: " + str(words))
                            raise ValueError(f"Required word {self.wordnum} does not exist")

                        value = 0.0
                        try:
                            value = int(words[self.wordnum])
                        except Exception:
                            self.annotations.append("Python error")
                            self.SetLastResultSummary(
                                f'IntCompare key="{self.key}" expected={self.gold} '
                                f'result=FAILED reason=parse_error')
                            if verbose:
                                warnings.warn("    Check failed :\n"
                                              + "    Check = " + self.__str__() + "\n"
                                              + "    line = " + line
                                              + "    words = " + str(words) + "\n"
                                              + "    info = Failed to convert word "
                                              + str(self.wordnum) + " to int.")
                            return False

                        if value == self.gold:
                            self.SetLastResultSummary(
                                f'IntCompare key="{self.key}" expected={self.gold} '
                                f'actual={value} result=PASSED')
                            return True
                        elif verbose:
                            warnings.warn("Check failed : " + self.__str__() + "\n"
                                          + line + "\n"
                                          + str(words))
                        self.SetLastResultSummary(
                            f'IntCompare key="{self.key}" expected={self.gold} '
                            f'actual={value} result=FAILED')

                if verbose:
                    warnings.warn('Check failed : key, "' + self.key + '", not found, wordnum = '
                                  + str(self.wordnum))
                self.SetLastResultSummary(
                    f'IntCompare key="{self.key}" expected={self.gold} '
                    f'result=FAILED reason=key_not_found')

        except FileNotFoundError as e:
            print(str(e))
            self.SetLastResultSummary(
                f'IntCompare key="{self.key}" expected={self.gold} '
                f'result=FAILED reason=file_not_found')
        except Exception as e:
            self.annotations.append("Python error")
            self.SetLastResultSummary(
                f'IntCompare key="{self.key}" expected={self.gold} '
                f'result=FAILED reason=python_error')
            if verbose:
                print(str(e))

        return False


class ErrorCodeCheck(Check):
    """Purely compares the error_code of the simulation. Useful for seeing if
      the program just runs (not crashing) and for checking that programs fail
      in a sane way"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.error_code: int = 0

        CheckUnknownFields(params, {"type", "error_code"}, message_prefix)
        RequireInteger(params, "error_code", message_prefix)

        self.error_code = params["error_code"]

    def __str__(self):
        return f'error_code="{self.error_code}"'

    def PerformCheck(self, filename, errorcode, verbose: bool):
        try:
            if errorcode == self.error_code:
                self.SetLastResultSummary(
                    f'ErrorCode expected={self.error_code} actual={errorcode} result=PASSED')
                return True

            if verbose:
                warnings.warn(f'Check failed: error_code, {self.error_code} vs {errorcode}')
            self.SetLastResultSummary(
                f'ErrorCode expected={self.error_code} actual={errorcode} result=FAILED')

        except Exception as e:
            self.annotations.append("Python error")
            self.SetLastResultSummary(
                f'ErrorCode expected={self.error_code} actual={errorcode} '
                f'result=FAILED reason=python_error')
            if verbose:
                warnings.warn(str(e))

        return False


# ===================================================================
class GoldFileCheck(Check):
    """Compares the output of a test against a gold-file"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.scope_keyword: str = ""
        self.candidate_filename: str = ""
        self.skiplines_top: int = 0
        self.check_numlines: int = 0

        CheckUnknownFields(params,
                           {"type", "scope_keyword", "candidate_filename",
                            "skiplines_top", "check_numlines"},
                           message_prefix)
        if "scope_keyword" in params:
            if not isinstance(params["scope_keyword"], str):
                warnings.warn(message_prefix + '"scope_keyword" field must be a string')
                raise ValueError
            self.scope_keyword = params["scope_keyword"]
        if "candidate_filename" in params:
            if not isinstance(params["candidate_filename"], str):
                warnings.warn(message_prefix + '"candidate_filename" field must be a string')
                raise ValueError
            self.candidate_filename = params["candidate_filename"]
        if "skiplines_top" in params:
            RequireInteger(params, "skiplines_top", message_prefix)
            self.skiplines_top = params["skiplines_top"]
        if "check_numlines" in params:
            RequireInteger(params, "check_numlines", message_prefix)
            self.check_numlines = params["check_numlines"]

    def __str__(self):
        if self.scope_keyword != "":
            return f'scope_begin_keyword="{self.scope_keyword}_BEGIN" ' + \
                f'scope_end_keyword="{self.scope_keyword}_END" '
        return "pure_diff"

    def PerformCheck(self, filename, errorcode, verbose: bool):
        self.SetLastResultSummary('GoldFile result=FAILED reason=unknown')
        try:
            outfiledir = pathlib.Path(os.path.dirname(filename) + "/")
            if self.candidate_filename != "":
                filename = str(outfiledir.parent.absolute()) + "/" + self.candidate_filename
            golddir = str(outfiledir.parent.absolute()) + "/gold/"

            goldfilename = os.path.splitext(os.path.basename(filename))[0] + ".gold"
            if self.candidate_filename != "":
                goldfilename = self.candidate_filename + ".gold"

            if not os.path.isfile(golddir + goldfilename):
                if verbose:
                    warnings.warn(f'Gold file {goldfilename} does not exist at \n'
                                  + f'{golddir + goldfilename}')
                self.annotations.append("Gold file missing")
                self.SetLastResultSummary(
                    f'GoldFile file="{goldfilename}" result=FAILED reason=gold_missing')
                return False

            lines_a, lines_b = self.GetFileLines(filename, golddir + goldfilename)

            if len(lines_a) == 0:
                if verbose:
                    print(f"no lines to compare in {filename}. "
                          + f"Maybe {self.scope_keyword}_BEGIN/_END was not found?")
                self.SetLastResultSummary(
                    f'GoldFile file="{goldfilename}" result=FAILED reason=no_candidate_lines')
                return False

            if len(lines_b) == 0:
                if verbose:
                    print(f"no lines to compare in {golddir + goldfilename}. "
                          + f"Maybe {self.scope_keyword}_BEGIN/_END was not found?")
                self.SetLastResultSummary(
                    f'GoldFile file="{goldfilename}" result=FAILED reason=no_gold_lines')
                return False

            if self.check_numlines == 0:
                diff = list(difflib.unified_diff(lines_a[self.skiplines_top:],
                                                 lines_b[self.skiplines_top:],
                                                 fromfile=filename,
                                                 tofile=golddir + goldfilename,
                                                 n=0  # Removes context
                                                 ))
            else:
                diff = list(difflib.unified_diff(
                    lines_a[self.skiplines_top:self.skiplines_top + self.check_numlines],
                    lines_b[self.skiplines_top:self.skiplines_top + self.check_numlines],
                    fromfile=filename,
                    tofile=golddir + goldfilename,
                    n=0  # Removes context
                ))

            if len(diff) == 0:
                self.SetLastResultSummary(
                    f'GoldFile file="{goldfilename}" diff_lines=0 result=PASSED')
                return True
            elif verbose:
                print("diff-begin")
                for line in diff:
                    print(line, end='')
                print("diff-end")
            self.SetLastResultSummary(
                f'GoldFile file="{goldfilename}" diff_lines={len(diff)} result=FAILED')

        except Exception as e:
            self.annotations.append("Python error")
            self.SetLastResultSummary(
                f'GoldFile file="{filename}" result=FAILED reason=python_error')
            if verbose:
                warnings.warn(str(e))

        return False

    def GetFileLines(self, filename_a: str, filename_b: str):
        def GetRawLines(filename):
            with open(filename, "r", encoding='utf-8') as file:
                lines = file.readlines()
            return lines

        def ScopeFilterLines(input_lines: list, scope_keyword: str):
            lines = []
            read_gate_open = False
            for line in input_lines:
                if line.find(scope_keyword + "_BEGIN") >= 0:
                    read_gate_open = True

                if line.find(scope_keyword + "_END") >= 0:
                    read_gate_open = False

                if read_gate_open:
                    lines.append(line)

            return lines

        lines_a = GetRawLines(filename_a)
        lines_b = GetRawLines(filename_b)

        if self.scope_keyword != "":
            lines_a = ScopeFilterLines(lines_a, self.scope_keyword)
            lines_b = ScopeFilterLines(lines_b, self.scope_keyword)

        return lines_a, lines_b
