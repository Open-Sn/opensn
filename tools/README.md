# Contents 

This `tools` directory contains various tools for use by both developers and users.

Files located at the root of `tools`:
- `configure_dependencies.py`: Script for installing OpenSn dependencies
- [lua-input-style.md](./lua-input-style.md): Notes on how to format `lua` input files. 

### Tools for Developer

Files are located in `tools/developer`.

A description of the files or sub-directory is given in the table below.

| File or Directory  | Description                                                                                    |
|--------------------|------------------------------------------------------------------------------------------------|
| `lsan.supp`        | Leak sanitizer suppression file (used in git workflows)                                        |
| `ignorelist.txt`   | List of source-level entities that should be ignored by sanitizers (used in CMakepresets.json) |

### Scaling Studies

Files are located in `scaling`.

A description of the files or sub-directory is given in the table below.

| File or Directory | Description |
|-------------------| ----------- |
|  `scaling`        | OpenSn scaling studies (see `README` for each study for details) |

### Cross-Section Generation

Files are located in `xs`.

A description of the files or sub-directory is given in the table below.

| File or Directory | Description |
| ----------------- | ----------- |
| [XS README](./MGXS_README.md) | Description of the contents of the `xs` directory |



[Return to OpenSn](../README.md)
