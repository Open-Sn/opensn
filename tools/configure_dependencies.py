"""
This is a utility script for the download and installation of OpenSn
dependencies.

NOTES:
- If you see errors in building fblaslapack, try using a  different BLAS/LAPACK. For example, on
  MacOSX, the Apple-supplied ones are known to work (so you would just remove the option
  --download-fblaslapack=1)
- If you have set CC, CXX and FC in your environment, make sure they are CC=mpicc, CXX=mpicxx, and
  FC=mpifort.

NOTES for building on MacOS:
- you may need to `export MACOSX_DEPLOYMENT_TARGET=10.15` in your environment
"""

import os
import sys
import argparse
import shutil
import subprocess
import textwrap
from typing import Optional, TextIO

if sys.version_info.major < 3:
    raise Exception("Python version detected to be " + str(sys.version_info.major) + "."
                    + str(sys.version_info.minor) + " but is required to >= 3")


class TextColors:
    RED = "\033[31m"
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


error_beg = TextColors.RED
error_end = TextColors.ENDC

# Process commandline arguments
arguments_help = textwrap.dedent('''\
Run the dependency builder.
''')

parser = argparse.ArgumentParser(
    description="This is a utility script for the download and installation of OpenSn "
                + "dependencies.",
    epilog=arguments_help
)

parser.add_argument(
    "-d", "--directory",
    type=str, required=True, default=None,
    help="Directory into which the dependencies will be installed"
)

parser.add_argument(
    "-j", "--jobs",
    type=int, required=False, default=4,
    help="Allow N compile jobs at once"
)

parser.add_argument(
    "-v", "--verbose",
    type=bool, required=False, default=False,
    help="Controls verbose failure"
)

parser.add_argument(
    "-dl", "--download_only",
    type=bool, required=False, default=False,
    help="Controls verbose failure"
)

# Define versions and download paths
VERSION = 0
URL = 1
package_info = {
    "boost": [
        "1_86_0",
        "https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz"
    ],
    "lua": [
        "5.4.6",
        "https://www.lua.org/ftp/lua-5.4.6.tar.gz"
    ],
    "hdf5": [
        "1.14.3",
        "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/src/\
CMake-hdf5-1.14.3.tar.gz"
    ],
    "petsc": [
        "3.20.6",
        "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.20.6.tar.gz"
    ],
    "vtk": [
        "9.3.0",
        "https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz"
    ],
    "caliper": [
        "2.10.0",
        "https://github.com/LLNL/Caliper/archive/refs/tags/v2.10.0.tar.gz"
    ]
}


# Runs a subprocess
def ExecSub(command: str,
            out_log: Optional[TextIO] = subprocess.PIPE,
            err_log: Optional[TextIO] = subprocess.PIPE,
            env_vars=None):
    success = True

    if argv.verbose:
        print(f"command: {command}")

    result = subprocess.Popen(command,
                              stdout=out_log,
                              stderr=err_log,
                              shell=True,
                              env=env_vars)
    output, error = result.communicate()

    if result.returncode != 0:
        success = False

    outputstr = ""
    if output is not None:
        outputstr = output.decode('ascii')

    return success, error.decode("utf8"), outputstr


# Checks whether an executable exists
def CheckExecutableExists(thingname: str, thing: str):
    log_file.write(f"Checking for {thingname}:\n")
    log_file.flush()
    success, errorc, outstr = ExecSub(f"{thing} --version", out_log=log_file)
    if not success:
        log_file.write(errorc)
        raise RuntimeError(f"Valid {thingname} not found")
    else:
        log_file.write(f"{thingname} found\n\n")


# Downloads a package using curl
def DownloadPackage(url, pkg, ver):
    pkgdir = f"{install_dir}/downloads"
    log_file.write(f"Downloading {pkg.upper()} {ver} to \"{pkgdir}\" "
                   + f"with command: curl {url}")
    download_cmd = "curl"
    output_cmd = f"-L --output {pkg}-{ver}.tar.gz"

    pkgdir_relative = os.path.relpath(pkgdir)

    print(f"Downloading {pkg.upper()} {ver} to \"{pkgdir_relative}\" "
          + f"with command:\n{download_cmd} {url} {output_cmd}", end="", flush=True)

    already_there = False
    if not os.path.exists(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz"):
        ExecSub(f"{download_cmd} {url} {output_cmd}", out_log=log_file)

        if os.path.exists(f"{install_dir}/downloads/{pkg.upper()}-{ver}.tar.gz"):
            item = f"{pkg.upper()}-{ver}.tar.gz"
            shutil.move(f"{item}", f"{item.lower()}")
    else:
        already_there = True

    if os.path.exists(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz"):
        there = "Already downloaded" if already_there else ""
        print(f" {TextColors.OKGREEN}Success{TextColors.ENDC} "
              + f"{TextColors.OKCYAN}{there}{TextColors.ENDC}")
        log_file.write(f" Success {there}\n")
        return True
    else:
        print(error_beg + " Failure" + error_end)
        log_file.write(" Failure\n")
        return False


# Makes a directory if it does not exist yet
def MakeDirectory(dirpath: str):
    if not os.path.isdir(dirpath):
        os.mkdir(dirpath)


# Makes the necessary system calls to extract a tar.gz file
def ExtractPackage(pkg, ver):
    print(f"Extracting package with command tar -zxf {pkg}-{ver}.tar.gz")

    success, err, outstr = ExecSub(f"tar -zxf {pkg}-{ver}.tar.gz", log_file)

    if not success:
        raise RuntimeError(err)


def InstallBoostPackage(pkg: str,
                        ver: str,
                        gold_file: str):
    package_log_filename = f"{install_dir}/logs/{pkg}_log.txt"
    pkg_install_dir = f"{install_dir}"

    # Copy the downloaded tarball to the source directory
    shutil.copy(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz",
                f"{install_dir}/src/{pkg}-{ver}.tar.gz")

    os.chdir(f"{install_dir}/src")

    # Check if Boost headers are already installed
    if not os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        ExtractPackage(pkg, ver)

        package_log_file = open(package_log_filename, "w")

        print(f"Installing {pkg.upper()} {ver} headers in \"{pkg_install_dir}/include/boost\"",
              flush=True)
        log_file.write(f"Installing {pkg.upper()} {ver} headers to\
                       \"{pkg_install_dir}/include/boost\"")
        log_file.write(f" See {package_log_filename}\n")
        log_file.flush()

        os.chdir(f"{pkg}_{ver}")

        # Set up the Boost include directory
        boost_include_dir = os.path.join(pkg_install_dir, "include")
        os.makedirs(boost_include_dir, exist_ok=True)

        # Copy header files to the Boost include directory
        command = f"cp -r boost {boost_include_dir}"
        success, err, outstr = ExecSub(command, out_log=package_log_file)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to install {pkg} headers")

        package_log_file.close()
    else:
        print(f"{pkg} already installed")

    os.chdir(install_dir)

    if os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        return True
    else:
        return False


# Install Lua
def InstallLuaPackage(pkg: str,
                      ver: str,
                      gold_file: str):
    package_log_filename = f"{install_dir}/logs/{pkg}_log.txt"
    pkg_install_dir = f"{install_dir}"

    shutil.copy(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz",
                f"{install_dir}/src/{pkg}-{ver}.tar.gz")

    os.chdir(f"{install_dir}/src")

    # Check if it is installed already
    if not os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        ExtractPackage(pkg, ver)

        package_log_file = open(package_log_filename, "w")

        print(f"Configuring {pkg.upper()} {ver} in \"{os.getcwd()}\"", flush=True)
        log_file.write(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"")
        log_file.write(f" See {package_log_filename}\n")
        log_file.flush()

        env_vars = os.environ.copy()
        if len(os.getenv("CC")) == 0:
            env_vars["CC"] = shutil.which("mpicc")
        if len(os.getenv("CXX")) == 0:
            env_vars["CXX"] = shutil.which("mpicxx")
        if len(os.getenv("FC")) == 0:
            env_vars["FC"] = shutil.which("mpifort")

        os_tag = "linux"
        if "Darwin" in os.uname():
            os_tag = "macosx"

        os.chdir(f"{pkg}-{ver}")

        command = f"make {os_tag} MYCFLAGS=-fPIC MYLIBS=-lncurses -j{argv.jobs}"
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to build {pkg}")

        command = f"make install INSTALL_TOP={pkg_install_dir}"
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to install {pkg}")

        package_log_file.close()
    else:
        print(f"{pkg} already installed")

    os.chdir(install_dir)

    if os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        return True
    else:
        return False


# Install PETSc
def InstallPETSc(pkg: str, ver: str, gold_file: str):
    package_log_filename = f"{install_dir}/logs/{pkg}_log.txt"
    pkg_install_dir = f"{install_dir}"

    shutil.copy(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz",
                f"{install_dir}/src/{pkg}-{ver}.tar.gz")

    os.chdir(f"{install_dir}/src")

    # Check if it is installed already
    if not os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        ExtractPackage(pkg, ver)

        package_log_file = open(package_log_filename, "w")

        print(f"Configuring {pkg.upper()} {ver} in \"{os.getcwd()}\"", flush=True)
        log_file.write(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"")
        log_file.write(f" See {package_log_filename}\n")
        log_file.flush()

        env_vars = os.environ.copy()
        if len(os.getenv("CC")) == 0:
            env_vars["CC"] = shutil.which("mpicc")
        if len(os.getenv("CXX")) == 0:
            env_vars["CXX"] = shutil.which("mpicxx")
        if len(os.getenv("FC")) == 0:
            env_vars["FC"] = shutil.which("mpifort")

        os.chdir(f"{pkg}-{ver}")
        command = f"""./configure --prefix={pkg_install_dir} \\
--download-hypre=1 \\
--with-ssl=0 \\
--with-debugging=0 \\
--with-pic=1 \\
--with-shared-libraries=1 \\
--download-bison=1 \\
--download-fblaslapack=1 \\
--download-metis=1 \\
--download-parmetis=1 \\
--download-superlu_dist=1 \\
--download-ptscotch=1 \\
--with-cxx-dialect=C++11 \\
--with-64-bit-indices \\
CC=$CC CXX=$CXX FC=$FC \\
FFLAGS="-fallow-argument-mismatch" \\
FCFLAGS="-fallow-argument-mismatch" \\
F90FLAGS="-fallow-argument-mismatch" \\
F77FLAGS="-fallow-argument-mismatch" \\
COPTFLAGS="-O3" \\
CXXOPTFLAGS="-O3" \\
FOPTFLAGS="-O3 -fallow-argument-mismatch" \\
PETSC_DIR={install_dir}/src/{pkg}-{ver}"""
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(
                f"Failed to configure {pkg}. See {package_log_filename} for details.")

        command = f"{make_command} all -j{argv.jobs}"
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to build {pkg}")

        command = f"{make_command} install"
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to install {pkg}")

        package_log_file.close()
    else:
        print(f"{pkg} already installed")

    os.chdir(install_dir)

    if os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        return True
    else:
        return False


# Install VTK
def InstallVTK(pkg: str, ver: str, gold_file: str):
    package_log_filename = f"{install_dir}/logs/{pkg}_log.txt"
    pkg_install_dir = f"{install_dir}"

    shutil.copy(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz",
                f"{install_dir}/src/{pkg}-{ver}.tar.gz")

    os.chdir(f"{install_dir}/src")

    # Check if it is installed already
    if not os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        ExtractPackage(pkg, ver)

        package_log_file = open(package_log_filename, "w")

        print(f"Configuring {pkg.upper()} {ver} in \"{os.getcwd()}\"", flush=True)
        log_file.write(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"")
        log_file.write(f" See {package_log_filename}\n")
        log_file.flush()

        env_vars = os.environ.copy()
        if len(os.getenv("CC")) == 0:
            env_vars["CC"] = shutil.which("mpicc")
        if len(os.getenv("CXX")) == 0:
            env_vars["CXX"] = shutil.which("mpicxx")
        if len(os.getenv("FC")) == 0:
            env_vars["FC"] = shutil.which("mpifort")

        build_dir = f"{install_dir}/src/{pkg.upper()}-{ver}/build"
        MakeDirectory(build_dir)
        os.chdir(build_dir)

        command = f""" cmake -DCMAKE_INSTALL_PREFIX={pkg_install_dir} \\
-DBUILD_SHARED_LIBS=ON \\
-DVTK_USE_MPI=ON \\
-DVTK_GROUP_ENABLE_StandAlone=WANT \\
-DVTK_GROUP_ENABLE_Rendering=DONT_WANT \\
-DVTK_GROUP_ENABLE_Imaging=DONT_WANT \\
-DVTK_GROUP_ENABLE_Web=DONT_WANT \\
-DVTK_GROUP_ENABLE_Qt=DONT_WANT \\
-DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON \\
-DCMAKE_BUILD_TYPE=Release \\
../
"""
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to configure {pkg}")

        command = f"{make_command} -j{argv.jobs}"
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to build {pkg}")

        command = f"{make_command} install"
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to install {pkg}")

        package_log_file.close()
    else:
        print(f"{pkg} already installed")

    os.chdir(install_dir)

    if os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        return True
    else:
        return False


# Install Caliper
def InstallCaliper(pkg: str, ver: str, gold_file: str):
    package_log_filename = f"{install_dir}/logs/{pkg}_log.txt"
    pkg_install_dir = f"{install_dir}"

    shutil.copy(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz",
                f"{install_dir}/src/caliper-{ver}.tar.gz")

    os.chdir(f"{install_dir}/src")

    # Check if it is installed already
    if not os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        ExtractPackage(pkg, ver)

        package_log_file = open(package_log_filename, "w")

        print(f"Configuring {pkg.upper()} {ver} in \"{os.getcwd()}\"", flush=True)
        log_file.write(f"Configuring {pkg.upper()} {ver} in \"{os.getcwd()}\"")
        log_file.write(f" See {package_log_filename}\n")
        log_file.flush()

        env_vars = os.environ.copy()
        if len(os.getenv("CC")) == 0:
            env_vars["CC"] = shutil.which("mpicc")
        if len(os.getenv("CXX")) == 0:
            env_vars["CXX"] = shutil.which("mpicxx")
        if len(os.getenv("FC")) == 0:
            env_vars["FC"] = shutil.which("mpifort")

        build_dir = f"{install_dir}/src/Caliper-{ver}/build"
        MakeDirectory(build_dir)
        os.chdir(build_dir)

        command = f""" cmake -DCMAKE_INSTALL_PREFIX={pkg_install_dir} \\
-DWITH_MPI=ON \
-DWITH_KOKKOS=OFF \
../
"""
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to configure {pkg}")

        command = f"{make_command} -j{argv.jobs}"
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to build {pkg}")

        command = f"{make_command} install"
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to install {pkg}")

        package_log_file.close()
    else:
        print(f"{pkg} already installed")

    os.chdir(install_dir)

    if os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        return True
    else:
        return False


# Install HDF5
def InstallHDF5(pkg: str, ver: str, gold_file: str):
    package_log_filename = f"{install_dir}/logs/{pkg}_log.txt"
    pkg_install_dir = f"{install_dir}"

    shutil.copy(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz",
                f"{install_dir}/src/{pkg}-{ver}.tar.gz")

    os.chdir(f"{install_dir}/src")

    # Check if it is installed already
    if not os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        ExtractPackage(pkg, ver)

        package_log_file = open(package_log_filename, "w")

        print(f"Configuring {pkg.upper()} {ver} in \"{os.getcwd()}\"", flush=True)
        log_file.write(f"Configuring {pkg.upper()} {ver} in \"{os.getcwd()}\"")
        log_file.write(f" See {package_log_filename}\n")
        log_file.flush()

        env_vars = os.environ.copy()
        if len(os.getenv("CC")) == 0:
            env_vars["CC"] = shutil.which("mpicc")
        if len(os.getenv("CXX")) == 0:
            env_vars["CXX"] = shutil.which("mpicxx")
        if len(os.getenv("FC")) == 0:
            env_vars["FC"] = shutil.which("mpifort")

        build_dir = f"{install_dir}/src/CMake-{pkg}-{ver}/{pkg}-{ver}/build"
        MakeDirectory(build_dir)
        os.chdir(build_dir)

        command = f""" cmake -DCMAKE_INSTALL_PREFIX={pkg_install_dir} .."""
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to configure {pkg}")

        command = f"{make_command} -j{argv.jobs}"
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to build {pkg}")

        command = f"{make_command} install"
        success, err, outstr = ExecSub(
            command, out_log=package_log_file, env_vars=env_vars
        )
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")
            raise RuntimeError(f"Failed to install {pkg}")

        package_log_file.close()
    else:
        print(f"{pkg} already installed")

    os.chdir(install_dir)

    if os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        return True
    else:
        return False


try:
    argv = parser.parse_args()  # argv = argument values

    # Redefine make command for gnu
    make_command = "make"
    s, e, outstr = ExecSub(make_command + " --version")
    if outstr.find("GNU Make") >= 0:
        make_command = "make OMAKE_PRINTDIR=make"

    # Check the test directory exists
    if not os.path.isdir(argv.directory):
        raise NotADirectoryError(argv.directory)

    install_dir = os.path.abspath(argv.directory)

    MakeDirectory(f"{install_dir}/src")
    MakeDirectory(f"{install_dir}/logs")

    print(f"Installation directory identified: {install_dir}")
    print()
    log_file = open(f"{install_dir}/logs/log.txt", 'w')
    log_file.write("***** OpenSn dependency configurator script *****\n\n")

    log_file.write("Packages in the dependency list:\n")
    for pkg in package_info:
        log_file.write(f"{pkg:9s}= {package_info[pkg][VERSION]}\n")
    log_file.write("\n")

    log_file.write("Download paths:\n")
    for pkg in package_info:
        log_file.write(f"{pkg:9s}= {package_info[pkg][URL]}\n")
    log_file.write("\n")

    # Check system environment

    # Check for download utils
    log_file.write("Checking for curl: ")
    if shutil.which("curl") is None:
        log_file.write("curl not found. ")
        raise Exception("Missing Curl utility. Unable to download dependencies.")
    else:
        log_file.write("curl found\n")

    # Check mpicc, mpicxx, mpifort
    mpicc = os.getenv("CC") if len(os.getenv("CC")) > 0 else "mpicc"
    mpicxx = os.getenv("CXX") if len(os.getenv("CXX")) > 0 else "mpicxx"
    mpifort = os.getenv("FC") if len(os.getenv("FC")) > 0 else "mpifort"
    CheckExecutableExists("mpicc", mpicc)
    CheckExecutableExists("mpicxx", mpicxx)
    CheckExecutableExists("mpifort", mpifort)

    # Check cmake exists
    CheckExecutableExists("cmake", "cmake")

    # Creating directories
    MakeDirectory(f"{install_dir}/downloads")

    # Download packages
    os.chdir(f"{install_dir}/downloads")
    dl_errors = []
    for pkg in package_info:
        if not DownloadPackage(package_info[pkg][URL], pkg,
                               package_info[pkg][VERSION]):
            dl_errors.append(pkg)
    os.chdir(f"{install_dir}")

    if len(dl_errors) > 0:
        log_file.write(
            "\nFailed to download the following packages:\n")
        for pkg in dl_errors:
            log_file.write(pkg + " ")
        log_file.write("\n")
        print(error_beg + f"Failed to download {len(dl_errors)} package(s)\n"
              + "Check that the url exists or that the platform you are on actually allows it. You "
              + "could also try to manually download these packages from the URLs below:\n")
        for pkg in dl_errors:
            print(f"{pkg:9s}= {package_info[pkg][URL]}")
        print(error_end)
        exit(1)

    if argv.download_only:
        exit(0)

    print()

    #  Collective installation point
    log_file.write("\nInstalling packages:\n")
    log_file.flush()
    install_errors = []

    env_script_name = f"{install_dir}/bin/set_opensn_env.sh"

    for pkg in package_info:
        log_file.write(f"{pkg}\n")
        log_file.flush()
        ver = package_info[pkg][VERSION]
        success = False
        if pkg == 'boost':
            success = InstallBoostPackage(pkg, ver, gold_file="include/boost")
        elif pkg == 'lua':
            success = InstallLuaPackage(pkg, ver, gold_file="lib/liblua.a")
        elif pkg == 'petsc':
            success = InstallPETSc(pkg, ver, gold_file="include/petsc")
        elif pkg == 'vtk':
            major, minor, patch = ver.split('.')
            success = InstallVTK(pkg, ver, gold_file=f"include/vtk-{major}.{minor}")
        elif pkg == 'caliper':
            major, minor, patch = ver.split('.')
            success = InstallCaliper(pkg, ver, gold_file="include/caliper/cali.h")
        elif pkg == 'hdf5':
            major, minor, patch = ver.split('.')
            success = InstallHDF5(pkg, ver, gold_file="include/hdf5.h")
        else:
            print(f"No build rules for {pkg}")

        if not success:
            install_errors.append(pkg)

    if len(install_errors) > 0:
        log_file.write(
            "\nFailed to install the following packages:\n")
        for pkg in install_errors:
            log_file.write(pkg + " ")
        log_file.write("\n")
        print(error_beg + f"Failed to install {len(install_errors)} package(s)\n"
                          "Check the package's associated log file, i.e., "
                          "PACKAGE/package_log.txt for information as to why the install "
                          "failed. You could also try to manually install these packages")
        print(error_end)
        exit(1)

    # Create envvars file
    env_script_file = open(env_script_name, "w")

    petsc_dir = f"{install_dir}"
    vtk_dir = f"{install_dir}"

    env_script_file.write(
        f'export CMAKE_PREFIX_PATH="{install_dir}"'
        f'${{CMAKE_PREFIX_PATH:+:${{CMAKE_PREFIX_PATH}}}}\n')

    if os.path.exists(f"{vtk_dir}/lib64"):
        env_script_file.write(f'export LD_LIBRARY_PATH="{vtk_dir}/lib64":$LD_LIBRARY_PATH\n')
    else:
        env_script_file.write(f'export LD_LIBRARY_PATH="{vtk_dir}/lib":$LD_LIBRARY_PATH\n')

    env_script_file.close()

    ExecSub(f"chmod u+x {env_script_name}", log_file)
    log_file.close()

    print("\n########## OpenSn dependency install complete ##########")
    print("\nWhen opening OpenSn in an IDE, the following environment variables need to be set:\n")

    print(f'CMAKE_PREFIX_PATH="{install_dir}"${{CMAKE_PREFIX_PATH:+:${{CMAKE_PREFIX_PATH}}}}')
    print()
    print(TextColors.WARNING
          + "When compiling OpenSn, in a terminal, the following environment variables need to be "
          + "set:\n" + TextColors.ENDC)

    print(
        f'export CMAKE_PREFIX_PATH="{install_dir}"${{CMAKE_PREFIX_PATH:+:${{CMAKE_PREFIX_PATH}}}}')
    if os.path.exists(f"{vtk_dir}/lib64"):
        print(f'export LD_LIBRARY_PATH="{vtk_dir}/lib":$LD_LIBRARY_PATH')
    else:
        print(f'export LD_LIBRARY_PATH="{vtk_dir}/lib64":$LD_LIBRARY_PATH')

    print()
    print("To set these terminal environment variables automatically, execute:")
    print(f"    $ source {env_script_name}\n")

except Exception as e:
    print(f"{TextColors.RED}{e}{TextColors.ENDC}")
