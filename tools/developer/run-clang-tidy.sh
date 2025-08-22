#!/bin/bash
repo=$(git rev-parse --show-toplevel)
run-clang-tidy -p build                                              \
               -header-filter="^${repo}/(framework|modules|python)/" \
               "${repo}/modules" "${repo}/python"
