#!/usr/bin/env python

import subprocess

try:
    output = subprocess.check_output(["git", "clang-format", "--need_reformat"])
except subprocess.CalledProcessError:
    print("You need to run 'git clang-format' and stage the files again before commit")
    exit(1)
else:
    exit(0)
