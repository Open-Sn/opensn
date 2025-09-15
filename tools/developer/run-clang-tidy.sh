#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: run-clang-tidy.sh [-h] [FILENAME]

OpenSn-specific clang-tidy script.

Runs clang-tidy on the entire OpenSn repository or on a single source file.

This command should be run from the root of the OpenSn repository.

For example, to analyze the entire repository, run:
tools/developer/run-clang-tidy.sh

To analyze a single source (.cc) file, run:
tools/developer/run-clang-tidy.sh framework/utils/timer.cc:

Options:
  -h, --help    Show this help and exit.

EOF
}

# Parse command-line arguments
positional=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    -*) echo "Unknown option: $1" >&2; usage; exit 2 ;;
    *) positional+=( "$1" ); shift; continue ;;
  esac
done
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

repo=$(git rev-parse --show-toplevel)
header_filter="^${repo}/(framework|modules|python)/"

if (( ${#positional[@]} == 1 )); then
  file="${positional[0]}"
  if [[ -f build/compile_commands.json ]]; then
    if ! grep -Fq -- "$file" build/compile_commands.json; then
      echo "Error: '$file' not found in build/compile_commands.json." >&2
      echo "       Ensure it’s part of a target and rebuild the DB." >&2
      exit 1
    fi
  fi
  run-clang-tidy -p build -header-filter="$header_filter" "$file"
else
  run-clang-tidy -p build -header-filter="$header_filter" "${repo}/modules" "${repo}/python"
fi
