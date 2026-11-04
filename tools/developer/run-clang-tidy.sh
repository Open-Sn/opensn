#!/usr/bin/env bash

use_color=0
positional=()
before_opt=()
arg_opt=()
removed_opt=()
clang_tidy_args=()
had_issues=0
tmp_log="$(mktemp -t clang_tidy.XXXXXX)"
cleanup() { rm -f "$tmp_log"; }
trap cleanup EXIT

set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: run-clang-tidy.sh [-h] [-use-color] [-extra-arg-before=ARG] [-extra-arg=ARG] [-removed-arg=ARG] [FILENAME]

OpenSn-specific clang-tidy script. Requires LLVM 19 or newer.

Runs `run-clang-tidy` on every OpenSn source file in build/compile_commands.json: all files in
the repository except vendored code in external/ and generated files in the build directory. You
can specify a single source file instead of running on the entire repository.

Examples (run from repository root):
  Analyze everything:
    tools/developer/run-clang-tidy.sh

  Analyze one file:
    tools/developer/run-clang-tidy.sh framework/utils/timer.cc

Options:
  -h, --help               Show this help text and exit.
  -use-color, --use-color  Enable colored output when supported.
  -extra-arg-before=ARG    Forward ARG before the file list (passed verbatim to llvm's run-clang-tidy).
  -extra-arg=ARG           Forward ARG after the file list.
  -removed-arg=ARG         Remove ARG from the compile command (passed verbatim to llvm's
                           run-clang-tidy; needs LLVM 20 or newer).
EOF
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    -extra-arg-before=*)
      before_opt+=( -extra-arg-before="${1#*=}" ); shift; continue ;;
    -extra-arg-before)
      shift; before_opt+=( -extra-arg-before="${1:-}" ); [[ -n "${1:-}" ]] && shift; continue ;;
    -extra-arg=*)
      arg_opt+=( -extra-arg="${1#*=}" ); shift; continue ;;
    -extra-arg)
      shift; arg_opt+=( -extra-arg="${1:-}" ); [[ -n "${1:-}" ]] && shift; continue ;;
    -removed-arg=*)
      removed_opt+=( -removed-arg="${1#*=}" ); shift; continue ;;
    -removed-arg)
      shift; removed_opt+=( -removed-arg="${1:-}" ); [[ -n "${1:-}" ]] && shift; continue ;;
    -use-color|--use-color)
      use_color=1; shift; continue ;;
    --) shift; break ;;
    -*) echo "Unknown option: $1" >&2; usage; exit 2 ;;
    *) positional+=( "$1" ); shift; continue ;;
  esac
done

# Add color flag if requested
if (( use_color )); then
  clang_tidy_args+=( -use-color )
fi

# Collect anything after '--'
if [[ $# -gt 0 ]]; then positional+=( "$@" ); fi

if (( ${#positional[@]} > 1 )); then
  echo "Error: too many arguments" >&2
  usage
  exit 2
fi

if ! command -v run-clang-tidy >/dev/null 2>&1; then
  echo "Error: 'run-clang-tidy' not found." >&2
  echo "       Install clang-tidy (LLVM) or add 'run-clang-tidy' to your PATH." >&2
  exit 127
fi
# -exclude-header-filter, used below to skip vendored and generated headers, needs LLVM 19 or newer.
if ! run-clang-tidy --help 2>&1 | grep -- '-exclude-header-filter' >/dev/null; then
  echo "Error: 'run-clang-tidy' does not support -exclude-header-filter. LLVM 19 or newer is required." >&2
  exit 1
fi
# -removed-arg needs LLVM 20 or newer.
if [[ ${#removed_opt[@]} -gt 0 ]] && ! run-clang-tidy --help 2>&1 | grep -- '-removed-arg' >/dev/null; then
  echo "Error: 'run-clang-tidy' does not support -removed-arg. LLVM 20 or newer is required." >&2
  exit 1
fi

# The repository is the one containing this script. The git toplevel of the current directory can
# be a different checkout, for example when this tree is a copy inside another repository.
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)
if [[ ! -d build ]]; then
  echo "Error: 'build' directory not found. Run from the repository root after configuring." >&2
  exit 1
fi
build_dir=$(cd build && pwd -P)

# Every file of this repository except vendored code and generated files is checked, so new
# source directories are covered without editing this script. Paths are escaped for the regular
# expressions below.
regex_escape() { printf '%s' "$1" | sed 's/[][\.*^$+?(){}|]/\\&/g'; }
# CMake can record the checkout or the build directory through a symbolic link, and clang-tidy
# matches the paths of the compilation database without resolving links. The filters therefore
# also accept the path recorded in CMakeCache.txt when it leads to the same directory.
cache_value() { sed -n "s/^$1:INTERNAL=//p" build/CMakeCache.txt 2>/dev/null | head -n 1 || true; }
path_re() {
  local re
  re=$(regex_escape "$1")
  if [[ -n "$2" && "$2" != "$1" && -d "$2" && "$(cd "$2" && pwd -P)" == "$1" ]]; then
    re="(${re}|$(regex_escape "$2"))"
  fi
  printf '%s' "$re"
}
repo_re=$(path_re "$repo" "$(cache_value CMAKE_HOME_DIRECTORY)")
build_re=$(path_re "$build_dir" "$(cache_value CMAKE_CACHEFILE_DIR)")
header_filter="^${repo_re}/"
exclude_header_filter="^(${repo_re}/external/|${build_re}/)"
# run-clang-tidy matches this against compilation-database paths with Python's re module.
source_filter="^(?!${build_re}/)${repo_re}/(?!external/)"


if (( ${#positional[@]} == 1 )); then
  file="${positional[0]}"
  if [[ -f build/compile_commands.json ]]; then
    if ! grep -Fq -- "$file" build/compile_commands.json; then
      echo "Error: '$file' not found in build/compile_commands.json." >&2
      echo "       Ensure it’s part of a target and rebuild the DB." >&2
      exit 1
    fi
  fi
  run-clang-tidy ${clang_tidy_args[@]:-} -p build -header-filter="$header_filter" \
                 -exclude-header-filter="$exclude_header_filter" ${before_opt[@]+"${before_opt[@]}"} \
                 ${arg_opt[@]+"${arg_opt[@]}"} ${removed_opt[@]+"${removed_opt[@]}"} "$file" \
                 -warnings-as-errors='*' | tee -a "$tmp_log"
else
  run-clang-tidy ${clang_tidy_args[@]:-} -p build -header-filter="$header_filter" \
                 -exclude-header-filter="$exclude_header_filter" ${before_opt[@]+"${before_opt[@]}"} \
                 ${arg_opt[@]+"${arg_opt[@]}"} ${removed_opt[@]+"${removed_opt[@]}"} "$source_filter" \
                 -warnings-as-errors='*' | tee -a "$tmp_log"
fi

# A filter that matches none of the compilation-database paths, for example when the build directory
# was configured from another checkout, would otherwise pass without analyzing anything.
if grep -E '^Running clang-tidy (in [0-9]+ threads )?for 0 files ' "$tmp_log" >/dev/null; then
  echo "Error: no files of build/compile_commands.json matched ${repo}." >&2
  exit 1
fi

# run-clang-tidy does not generate an exit status based on errors or warnings
# generated by clang-tidy. We'll have to do this ourselves by parsing the
# clang-tidy output.
if grep -E '^[^:]+:[0-9]+:[0-9]+: (warning|error): ' "$tmp_log" >/dev/null; then
  had_issues=1
fi

if (( ${had_issues} )); then
  echo "clang-tidy reported issues (treating warnings as errors)." >&2
  exit 1
fi
