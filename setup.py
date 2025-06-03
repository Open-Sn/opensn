#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat, February 22

Copyright (c) 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
"""

import os
import re
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    """
    CMake extension.

    This extension takes a source directory containing a ``CMakeLists.txt``
    file as an input to process.
    """

    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


class CMakeBuilder(build_ext):
    """
    CMake builder

    Build Python extension from a ``CMakeExtension`` in the build directory.
    Built binaries are then be copied to either the source directory (inplace
    build mode) or the library directory (pip install mode).
    """

    def run(self) -> None:
        """Execute the build extension process. Each extension is built in parallel."""
        # ensure that CMake is installed.
        try:
            subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError("CMake must be installed to build the extensions.")

        # build extensions
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext: CMakeExtension) -> None:
        """Execute CMake and build extension."""
        # get full path of the extension
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)
        extdir = ext_fullpath.parent.resolve()

        # get debug/release configuration
        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            "-DOPENSN_WITH_LUA=OFF",
            "-DOPENSN_WITH_PYTHON_MODULE=ON"
        ]
        if sys.platform.startswith("win"):
            cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]

        # get generator from environment if precised
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")
        if cmake_generator:
            cmake_args += [f"-G {cmake_generator}"]
        # get extra CMake arguments
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [
                extra_args
                for extra_args in os.environ["CMAKE_ARGS"].split(" ")
                if extra_args
            ]

        # special treatment for Darwin
        if sys.platform.startswith("darwin"):
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # setup build arguments
        build_args = []
        if sys.platform.startswith("win"):
            build_args += ["--config", cfg]

        # set number of threads for parallel build
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            if hasattr(self, "parallel") and self.parallel:
                build_args += [f"-j{self.parallel}"]
            else:
                build_args += [f"-j{os.cpu_count()}"]
        # get extra CMake arguments
        if "BUILD_ARGS" in os.environ:
            build_args += [
                extra_args
                for extra_args in os.environ["BUILD_ARGS"].split(" ")
                if extra_args
            ]

        # create temporary build directory
        build_temp = Path(self.build_temp) / ext.name
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        # configure and build
        subprocess.run(
            ["cmake", ext.sourcedir, *cmake_args], cwd=build_temp, check=True
        )
        subprocess.run(
            ["cmake", "--build", ".", *build_args], cwd=build_temp, check=True
        )


if __name__ == "__main__":
    setup(
        name="pyopensn",
        author="The OpenSn Authors",
        description="Python extension for OpenSn",
        packages=["pyopensn"],
        ext_modules=[CMakeExtension("pyopensn.__init__")],
        cmdclass={"build_ext": CMakeBuilder},
        install_requires=["mpi4py", "numpy"]
    )
