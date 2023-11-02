import os.path
import warnings
import re  # regular expressions
import pathlib
import difflib


class Check:
    """Base class for a Check data structure"""

    def __init__(self):
        self.annotations = []

    def __str__(self):
        return "Check base class"

    def PerformCheck(self, filename, errorcode, verbose: bool):
        warnings.warn("Unimplemented base class routine")
        return False

    def GetAnnotations(self):
        return self.annotations


# ===================================================================
class KeyValuePairCheck(Check):
    """Given a string key, checks the floating point value of the word
       immediately following the key"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.key: str = ""
        self.goldvalue: float = 0.0
        self.tol: float = 1.0
        self.skip_lines_until: str = ""

        if "key" not in params:
            warnings.warn(message_prefix + 'Missing "key" field')
            raise ValueError
        if "goldvalue" not in params:
            warnings.warn(message_prefix + 'Missing "goldvalue" field')
            raise ValueError
        if "tol" not in params:
            warnings.warn(message_prefix + 'Missing "tol" field')
            raise ValueError

        self.key = params["key"]
        self.goldvalue = params["goldvalue"]
        self.tol = params["tol"]

        if "skip_lines_until" in params:
            val = params["skip_lines_until"]
            if not isinstance(val, str):
                raise ValueError("Expected 'skip_lines_until' value to be str")
            self.skip_lines_until = val

    def __str__(self):
        return 'type="KeyValuePair", ' + f'key="{self.key}", ' + \
            f'goldvalue={self.goldvalue}, ' + \
            f'tol={self.tol}'

    def PerformCheck(self, filename, errorcode, verbose: bool):
        try:
            file = open(filename, "r")

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
                    postkey = line[(key_pos + len(self.key)):].strip()

                    words = re.split(r'\s+|,+]', postkey)
                    if len(words) == 0:
                        raise ValueError("word split failure: " + line)

                    try:
                        value = float(words[0])
                    except Exception as e:
                        self.annotations.append("Python error")
                        if verbose:
                            warnings.warn('Failed to convert word "'
                                          + words[0] + '" to float\n' +
                                          'post key: ' + postkey + "\n" +
                                          'words' + words.__str__())

                        return False

                    if abs(value - self.goldvalue) <= self.tol:
                        file.close()
                        return True
                    elif verbose:
                        print("Check failed : " + self.__str__() + "\n" + line)

            if verbose:
                print('Check failed : key, "' + self.key + '",not found ')

            file.close()
        except FileNotFoundError as e:
            warnings.warn(str(e))
        except Exception as e:
            self.annotations.append("Python error")
            if verbose:
                warnings.warn(str(e))

        return False


# ===================================================================
class StrCompareCheck(Check):
    """Given a key to identify a line, compares the specific word number against
       a golden value. The word is expected to be a string"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.key: str = ""
        self.wordnum: int = -1
        self.gold: str = ""
        self.skip_lines_until: str = ""

        if "key" not in params:
            warnings.warn(message_prefix + 'Missing "key" field')
            raise ValueError

        if "wordnum" in params:
            if "gold" not in params:
                warnings.warn(message_prefix + 'Missing "gold" field')
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
        try:
            file = open(filename, "r")

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
                        file.close()
                        return True
                    words = re.split(r'\s+|,+|=+', line.rstrip())

                    if len(words) <= self.wordnum:
                        warnings.warn("word count: " + str(len(words)) +
                                      "\nline: " +
                                      line.rstrip() +
                                      "\nwords: " + words.__str__())
                        raise ValueError(f"Required word {self.wordnum} does not exist")

                    value = words[self.wordnum]

                    if value == self.gold:
                        file.close()
                        return True
                    elif verbose:
                        warnings.warn("Check failed : " + self.__str__() + "\n" +
                                      line + "\n" + words.__str__())

            if verbose:
                warnings.warn('Check failed : key, "' + self.key +
                              '",not found ' + str(self.wordnum))

            file.close()
        except FileNotFoundError as e:
            print(str(e))
        except Exception as e:
            self.annotations.append("Python error")
            if verbose:
                print(str(e))

        return False


# ===================================================================
class FloatCompareCheck(Check):
    """Given a key to identify a line, compares the specific word number against
       a golden value. The word is expected to be a float"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.key: str = ""
        self.wordnum: int = 0
        self.gold: float = 0.0
        self.tol: float = 1.0e-6
        self.skip_lines_until: str = ""

        if "key" not in params:
            warnings.warn(message_prefix + 'Missing "key" field')
            raise ValueError
        if "wordnum" not in params:
            warnings.warn(message_prefix + 'Missing "wordnum" field')
            raise ValueError
        if "gold" not in params:
            warnings.warn(message_prefix + 'Missing "gold" field')
            raise ValueError
        if "tol" not in params:
            warnings.warn(message_prefix + 'Missing "tol" field')
            raise ValueError

        self.key = params["key"]
        self.wordnum = params["wordnum"]
        self.gold = params["gold"]
        self.tol = params["tol"]

        if "skip_lines_until" in params:
            val = params["skip_lines_until"]
            if not isinstance(val, str):
                raise ValueError("Expected 'skip_lines_until' value to be str")
            self.skip_lines_until = val

    def __str__(self):
        return 'type="FloatCompare", ' + f'key="{self.key}", ' + \
            f'wordnum={self.wordnum}, ' + \
            f'gold={self.gold} '

    def PerformCheck(self, filename, errorcode, verbose: bool):
        try:
            file = open(filename, "r")

            skip_lines = False
            perform_skip_check = bool(self.skip_lines_until != "")
            if perform_skip_check:
                skip_lines = True

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
                    words = re.split(r'\s+|,+|=+', line.rstrip())

                    if len(words) <= self.wordnum:
                        warnings.warn("word count: " + str(len(words)) +
                                      "\nline: " +
                                      line + "\n" +
                                      "\nwords: " + words.__str__())
                        raise ValueError(f"Required word {self.wordnum} does not exist")

                    value = 0.0
                    try:
                        value = float(words[self.wordnum])
                    except Exception as e:
                        self.annotations.append("Python error")
                        warnings.warn("    Check failed :\n" +
                                      "    Check = " + self.__str__() + "\n" +
                                      "    line = " + line +
                                      "    words = " + words.__str__() + "\n" +
                                      "    info = Failed to convert word " +
                                      str(self.wordnum) + " to float.")
                        return False

                    if abs(value - self.gold) <= self.tol:
                        file.close()
                        return True
                    elif verbose:
                        warnings.warn("Check failed : " + self.__str__() + "\n" +
                                      line + "\n" + words.__str__())

            if verbose:
                warnings.warn('Check failed : key, "' + self.key +
                              '",not found ' + str(self.wordnum))

            file.close()
        except FileNotFoundError as e:
            print(str(e))
        except Exception as e:
            self.annotations.append("Python error")
            if verbose:
                print(str(e))


        return False


# ===================================================================
class IntCompareCheck(Check):
    """Given a key to identify a line, compares the specific word number against
       a golden value. The word is expected to be a int"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.key: str = ""
        self.wordnum: int = 0
        self.gold: int = 0
        self.skip_lines_until: str = ""

        if "key" not in params:
            warnings.warn(message_prefix + 'Missing "key" field')
            raise ValueError
        if "wordnum" not in params:
            warnings.warn(message_prefix + 'Missing "wordnum" field')
            raise ValueError
        if "gold" not in params:
            warnings.warn(message_prefix + 'Missing "gold" field')
            raise ValueError

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
        try:
            file = open(filename, "r")

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
                    words = re.split(r'\s+|,+|=+', line.rstrip())

                    if len(words) <= self.wordnum:
                        warnings.warn("word count: " + str(len(words)) +
                                      "\nline: " +
                                      line + "\n" +
                                      "\nwords: " + words.__str__())
                        raise ValueError(f"Required word {self.wordnum} does not exist")

                    value = 0.0
                    try:
                        value = int(words[self.wordnum])
                    except Exception as e:
                        self.annotations.append("Python error")
                        warnings.warn("    Check failed :\n" +
                                      "    Check = " + self.__str__() + "\n" +
                                      "    line = " + line +
                                      "    words = " + words.__str__() + "\n" +
                                      "    info = Failed to convert word " +
                                      str(self.wordnum) + " to int.")
                        return False

                    if value == self.gold:
                        file.close()
                        return True
                    elif verbose:
                        warnings.warn("Check failed : " + self.__str__() + "\n" +
                                      line + "\n" + words.__str__())

            if verbose:
                warnings.warn('Check failed : key, "' + self.key +
                              '",not found, wordnum = ' + str(self.wordnum))

            file.close()
        except FileNotFoundError as e:
            print(str(e))
        except Exception as e:
            self.annotations.append("Python error")
            if verbose:
                print(str(e))

        return False


# ===================================================================
class ErrorCodeCheck(Check):
    """Purely compares the error_code of the simulation. Useful for seeing if
      the program just runs (not crashing) and for checking that programs fail
      in a sane way"""

    def __init__(self, params: dict, message_prefix: str):
        super().__init__()
        self.error_code: int = 0

        if "error_code" not in params:
            warnings.warn(message_prefix + 'Missing "error_code" field')
            raise ValueError

        self.error_code = params["error_code"]

    def __str__(self):
        return f'error_code="{self.error_code}"'

    def PerformCheck(self, filename, errorcode, verbose: bool):
        try:
            if errorcode == self.error_code:
                return True

            if verbose:
                warnings.warn('Check failed : error_code, ' +
                              f'{self.error_code} vs {errorcode}')

        except Exception as e:
            self.annotations.append("Python error")
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

        if "scope_keyword" in params:
            self.scope_keyword = params["scope_keyword"]
        if "candidate_filename" in params:
            self.candidate_filename = params["candidate_filename"]
        if "skiplines_top" in params:
            self.skiplines_top = params["skiplines_top"]

    def __str__(self):
        if self.scope_keyword != "":
            return f'scope_begin_keyword="{self.scope_keyword}_BEGIN" ' + \
                f'scope_end_keyword="{self.scope_keyword}_END" '
        else:
            return "pure_diff"

    def PerformCheck(self, filename, errorcode, verbose: bool):
        try:
            outfiledir = pathlib.Path(os.path.dirname(filename) + "/")
            if self.candidate_filename != "":
                filename = str(outfiledir.parent.absolute()) + "/" + \
                           self.candidate_filename
            golddir = str(outfiledir.parent.absolute()) + "/gold/"

            goldfilename = os.path.splitext(os.path.basename(filename))[0] + ".gold"
            if self.candidate_filename != "":
                goldfilename = self.candidate_filename + ".gold"

            if not os.path.isfile(golddir + goldfilename):
                if verbose:
                    warnings.warn(f'Gold file {goldfilename} does not exist at \n' +
                                  f'{golddir + goldfilename}')
                self.annotations.append("Gold file missing")
                return False

            lines_a, lines_b = self.GetFileLines(filename, golddir + goldfilename)

            if len(lines_a) == 0:
                if verbose:
                    print(f"no lines to compare in {filename}. " +
                          f"Maybe {self.scope_keyword}_BEGIN/_END was not found?")
                return False

            if len(lines_b) == 0:
                if verbose:
                    print(f"no lines to compare in {golddir + goldfilename}. " +
                          f"Maybe {self.scope_keyword}_BEGIN/_END was not found?")
                return False

            diff = list(difflib.unified_diff(lines_a[self.skiplines_top:],
                                             lines_b[self.skiplines_top:],
                                             fromfile=filename,
                                             tofile=golddir + goldfilename,
                                             n=0  # Removes context
                                             ))
            if len(diff) == 0:
                return True
            elif verbose:
                print("diff-begin")
                for line in diff:
                    print(line, end='')
                print("diff-end")

        except Exception as e:
            self.annotations.append("Python error")
            if verbose:
                warnings.warn(str(e))

        return False

    def GetFileLines(self, filename_a: str, filename_b: str):
        def GetRawLines(filename):
            file = open(filename, "r")
            lines = file.readlines()
            file.close()
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
