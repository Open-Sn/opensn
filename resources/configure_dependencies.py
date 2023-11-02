"""
This is a utility script for the download and installation of chi-tech
dependencies.
"""
import os
import sys
import argparse
import textwrap
import subprocess
import shutil
from typing import Optional, TextIO

if sys.version_info.major < 3:
    raise Exception("Python version detected to be " +
                    str(sys.version_info.major) + "." +
                    str(sys.version_info.minor) +
                    " but is required to >= 3")


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

# ========================================== Process commandline arguments
arguments_help = textwrap.dedent('''\
Run the dependency builder. 
''')

parser = argparse.ArgumentParser(
    description="This is a utility script for the download and installation" +
                "of chi-tech dependencies.",
    epilog=arguments_help
)

parser.add_argument(
    "-d", "--directory", default=None, type=str, required=True,
    help="Directory into which the dependencies will be installed"
)

parser.add_argument(
    "-j", "--jobs", default=4, type=int, required=False,
    help="Allow N compile jobs at once"
)

parser.add_argument(
    "-v", "--verbose", default=False, type=bool, required=False,
    help="Controls verbose failure"
)

parser.add_argument(
    "-dl", "--download_only", default=False, type=bool, required=False,
    help="Controls verbose failure"
)

argv = parser.parse_args()  # argv = argument values

# =========================================== Check stuff
# Check the test directory exists
if not os.path.isdir(argv.directory):
    raise NotADirectoryError(argv.directory)

install_dir = os.path.abspath(argv.directory)
print(f"Installation directory identified: {install_dir}")
print()
log_file = open(f"{install_dir}/log.txt", 'w')
log_file.write("***** ChiTech dependency configurator script *****\n\n")

# ========================================== Define versions and download paths
VERSION = 0
URL = 1
package_info = {
    "readline": ["8.0",
                 "ftp://ftp.gnu.org/gnu/readline/readline-8.0.tar.gz"],
    "ncurses": ["6.1",
                "https://invisible-mirror.net/archives/ncurses/ncurses-6.1.tar.gz"],
    "lua": ["5.4.6",
            "https://www.lua.org/ftp/lua-5.4.6.tar.gz"],
    "petsc": ["3.17.0",
              "https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.17.0.tar.gz"],
    "vtk": ["9.3.0.rc1",
            "https://www.vtk.org/files/release/9.3/VTK-9.3.0.rc1.tar.gz"]
}

log_file.write("Packages in the dependency list:\n")
for pkg in package_info:
    log_file.write(f"{pkg:9s}= {package_info[pkg][VERSION]}\n")
log_file.write("\n")

log_file.write("Download paths:\n")
for pkg in package_info:
    log_file.write(f"{pkg:9s}= {package_info[pkg][URL]}\n")
log_file.write("\n")


# ========================================== Define utilities
# ######################################
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


# ######################################
# Redefine make command for gnu
make_command = "make"
s, e, outstr = ExecSub(make_command + " --version")
if outstr.find("GNU Make") >= 0:
    make_command = "make OMAKE_PRINTDIR=gmake"


# ######################################
# Checks whether an executable exists
def CheckExecutableExists(thingname: str, thing: str):
    log_file.write(f"Checking for {thingname}:\n")
    log_file.flush()
    success, errorc, outstr = ExecSub(f"{thing} --version", out_log=log_file)
    if not success:
        log_file.write(errorc)
        print(error_beg + f"Valid {thingname} not found" + error_end)
        exit(1)
    else:
        log_file.write(f"{thingname} found\n\n")


# ###################################### Get package
# Downloads a package using either wget or curl
def DownloadPackage(downloader, url, pkg, ver):
    pkgdir = f"{install_dir}/downloads"
    log_file.write(f"Downloading {pkg.upper()} {ver} to \"{pkgdir}\" " +
                   f"with command: {downloader} {url}")
    download_cmd = f"{downloader}" if downloader == "wget" else f"{downloader}"
    output_cmd = "" if downloader == "wget" else f"--output {pkg}-{ver}.tar.gz"

    pkgdir_relative = os.path.relpath(pkgdir)

    print(f"Downloading {pkg.upper()} {ver} to \"{pkgdir_relative}\" " +
          f"with command:\n{download_cmd} {url} {output_cmd}", end="",
          flush=True)

    already_there = False
    if not os.path.exists(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz"):
        ExecSub(f"{download_cmd} {url} {output_cmd}", out_log=log_file)

        if os.path.exists(f"{install_dir}/downloads/{pkg.upper()}-{ver}.tar.gz"):
            item = f"{pkg.upper()}-{ver}.tar.gz"
            shutil.move(f"{item}", f"{item.lower()}")
    else:
        already_there = True

    if os.path.exists(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz"):
        there = f"Already downloaded" if already_there else ""
        print(f" {TextColors.OKGREEN}Success{TextColors.ENDC} {TextColors.OKCYAN}{there}{TextColors.ENDC}")
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
        print(err)
        exit(1)

    # Check if folder defaulted to upper
    if os.path.exists(f"{pkg.upper()}-{ver}") and not os.path.exists(f"{pkg}-{ver}"):
        success, err, outstr = ExecSub(f"mv {pkg.upper()}-{ver}/ {pkg}-{ver}/", log_file)

        if not success:
            print(err)
            exit(1)

    os.chdir(f"{pkg}-{ver}")


# ========================================== Check system environment
# Check for download utils
log_file.write("Checking for wget/curl: ")
downloader = ""
if shutil.which("wget") is None:
    log_file.write("wget not found. ")
    if shutil.which("curl") is None:
        log_file.write("curl not found. ")
        raise Exception("Neither wget nor curl found on this system. " +
                        "The script then has no means to download the " +
                        "dependencies.")
    else:
        log_file.write("curl found\n")
        downloader = "curl"
else:
    log_file.write("wget found\n")
    downloader = "wget"
# downloader = "curl"

# Check mpicc, mpicxx, mpifort
mpicc = os.getenv("CC") if len(os.getenv("CC")) > 0 else "mpicc"
mpicxx = os.getenv("CXX") if len(os.getenv("CXX")) > 0 else "mpicxx"
mpifort = os.getenv("FC") if len(os.getenv("FC")) > 0 else "mpifort"
CheckExecutableExists("mpicc", mpicc)
CheckExecutableExists("mpicxx", mpicxx)
CheckExecutableExists("mpifort", mpifort)

# Check cmake exists
CheckExecutableExists("cmake", "cmake")

# ========================================== Creating directories
MakeDirectory(f"{install_dir}/downloads")

# ========================================== Download packages
os.chdir(f"{install_dir}/downloads")
dl_errors = []
for pkg in package_info:
    if not DownloadPackage(downloader, package_info[pkg][URL], pkg,
                           package_info[pkg][VERSION]):
        dl_errors.append(pkg)
os.chdir(f"{install_dir}")

if len(dl_errors) > 0:
    log_file.write(
        "\nFailed to download the following packages:\n")
    for pkg in dl_errors:
        log_file.write(pkg + " ")
    log_file.write("\n")
    print(error_beg +
          f"Failed to download {len(dl_errors)} package(s)\n" +
          "Check that the url exists or that the platform you are on actually " +
          "actually allows it. You could also try to manually download these " +
          "packages from the url's below:\n")
    for pkg in dl_errors:
        print(f"{pkg:9s}= {package_info[pkg][URL]}")
    print(error_end)
    exit(1)

if argv.download_only:
    exit(0)

print()


# ========================================== Install command for ncurses and
# readline
def InstallPackage(pkg: str, ver: str, gold_file: str):
    package_log_filename = f"{install_dir}/{pkg.upper()}/{pkg}_log.txt"
    MakeDirectory(f"{install_dir}/{pkg.upper()}")
    pkg_install_dir = f"{install_dir}/{pkg.upper()}/{pkg}-{ver}/chi_build"

    shutil.copy(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz",
                f"{install_dir}/{pkg.upper()}/{pkg}-{ver}.tar.gz")

    os.chdir(f"{install_dir}/{pkg.upper()}/")

    # Check if it is installed already
    if not os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        ExtractPackage(pkg, ver)

        package_log_file = open(package_log_filename, "w")

        print(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"", flush=True)
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

        command = f"./configure --prefix={pkg_install_dir}"
        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        command = f"make -j{argv.jobs}"
        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        command = "make install"
        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        package_log_file.close()
    else:
        print(f"{pkg.capitalize()} already installed")

    os.chdir(install_dir)

    if os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        return True
    else:
        return False


# ========================================== Install command for lua
def InstallLuaPackage(pkg: str, ver: str, gold_file: str,
                      readline_install: str,
                      ncurses_install: str):
    package_log_filename = f"{install_dir}/{pkg.upper()}/{pkg}_log.txt"
    MakeDirectory(f"{install_dir}/{pkg.upper()}")
    pkg_install_dir = f"{install_dir}/{pkg.upper()}/{pkg}-{ver}/install"

    shutil.copy(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz",
                f"{install_dir}/{pkg.upper()}/{pkg}-{ver}.tar.gz")

    os.chdir(f"{install_dir}/{pkg.upper()}/")

    # Check if it is installed already
    if not os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        ExtractPackage(pkg, ver)

        package_log_file = open(package_log_filename, "w")

        print(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"", flush=True)
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

        lib_path = env_vars.get("LIBRARY_PATH", "")
        lib_path = f"{lib_path}:{readline_install}/lib:" \
                   f"{ncurses_install}/lib"
        env_vars["LIBRARY_PATH"] = lib_path

        c_path = env_vars.get('CPATH', "")
        c_path = f"{c_path}:{readline_install}/include"
        env_vars["CPATH"] = c_path

        os_tag = "linux"
        if "Darwin" in os.uname():
            os_tag = "macosx"

        command = f"make {os_tag} MYCFLAGS=-fPIC MYLIBS=-lncurses -j{argv.jobs}"
        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        command = "make local"
        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        package_log_file.close()
    else:
        print(f"{pkg.capitalize()} already installed")

    os.chdir(install_dir)

    if os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        return True
    else:
        return False


# ========================================== Install command for PETSc
def InstallPETSC(pkg: str, ver: str, gold_file: str):
    package_log_filename = f"{install_dir}/{pkg.upper()}/{pkg}_log.txt"
    MakeDirectory(f"{install_dir}/{pkg.upper()}")
    pkg_install_dir = f"{install_dir}/{pkg.upper()}/{pkg}-{ver}-install"

    shutil.copy(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz",
                f"{install_dir}/{pkg.upper()}/{pkg}-{ver}.tar.gz")

    os.chdir(f"{install_dir}/{pkg.upper()}/")

    # Check if it is installed already
    if not os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        ExtractPackage(pkg, ver)

        package_log_file = open(package_log_filename, "w")

        print(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"", flush=True)
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

        command = f"./configure --prefix={pkg_install_dir}"
        command += """ \\
--download-hypre=1  \\
--with-ssl=0  \\
--with-debugging=0  \\
--with-pic=1  \\
--with-shared-libraries=1  \\
--download-fblaslapack=1  \\
--download-metis=1  \\
--download-parmetis=1  \\
--download-superlu_dist=1  \\
--with-cxx-dialect=C++11  \\
--with-64-bit-indices  \\
CC=$CC CXX=$CXX FC=$FC  \\
CFLAGS='-fPIC -fopenmp'  \\
CXXFLAGS='-fPIC -fopenmp'  \\
FFLAGS='-fPIC -fopenmp'  \\
FCFLAGS='-fPIC -fopenmp'  \\
F90FLAGS='-fPIC -fopenmp'  \\
F77FLAGS='-fPIC -fopenmp'  \\
COPTFLAGS='-O3 -march=native -mtune=native'  \\
CXXOPTFLAGS='-O3 -march=native -mtune=native'  \\
FOPTFLAGS='-O3 -march=native -mtune=native'  \\"""
        command += f"PETSC_DIR={install_dir}/{pkg.upper()}/{pkg}-{ver}"

        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        command = f"{make_command} all -j{argv.jobs}"
        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        command = f"{make_command} install"
        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        package_log_file.close()
    else:
        print(f"{pkg.capitalize()} already installed")

    os.chdir(install_dir)

    if os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        return True
    else:
        return False


# ========================================== Install vtk
def InstallVTK(pkg: str, ver: str, gold_file: str):
    package_log_filename = f"{install_dir}/{pkg.upper()}/{pkg}_log.txt"
    MakeDirectory(f"{install_dir}/{pkg.upper()}")
    pkg_install_dir = f"{install_dir}/{pkg.upper()}/{pkg}-{ver}-install"

    shutil.copy(f"{install_dir}/downloads/{pkg}-{ver}.tar.gz",
                f"{install_dir}/{pkg.upper()}/{pkg}-{ver}.tar.gz")

    os.chdir(f"{install_dir}/{pkg.upper()}/")

    # Check if it is installed already
    if not os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        ExtractPackage(pkg, ver)

        package_log_file = open(package_log_filename, "w")
        MakeDirectory(f"{install_dir}/{pkg.upper()}/{pkg}-{ver}/chi_build")

        print(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"", flush=True)
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

        os.chdir("chi_build")

        command = f"cmake -DCMAKE_INSTALL_PREFIX={pkg_install_dir} "
        command += """ \\
-DBUILD_SHARED_LIBS:BOOL=ON \\
-DVTK_Group_MPI:BOOL=ON \\
-DVTK_GROUP_ENABLE_Qt=NO \\
-DVTK_GROUP_ENABLE_Rendering=NO \\
-DVTK_GROUP_ENABLE_Imaging=NO \\
-DVTK_GROUP_ENABLE_StandAlone=WANT \\
-DVTK_GROUP_ENABLE_Web=NO \\
-DVTK_BUILD_TESTING:BOOL=OFF \\
-DCMAKE_BUILD_TYPE=Release \\
-DCMAKE_CXX_FLAGS=-std=c++11 \\
../
"""
        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        command = f"{make_command} -j{argv.jobs}"
        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        command = f"{make_command} install"
        success, err, outstr = ExecSub(command, out_log=package_log_file, env_vars=env_vars)
        if not success:
            print(command, err)
            log_file.write(f"{command}\n{err}\n")
            package_log_file.write(f"{command}\n{err}\n")

        package_log_file.close()
    else:
        print(f"{pkg} already installed")

    os.chdir(install_dir)

    if os.path.exists(f"{pkg_install_dir}/{gold_file}"):
        return True
    else:
        return False


# ========================================== Collective install point
log_file.write("\nInstalling packages:\n")
log_file.flush()
install_errors = []
for pkg in package_info:
    log_file.write(f"{pkg}\n")
    log_file.flush()
    ver = package_info[pkg][VERSION]
    success = False
    if pkg == 'readline':
        success = InstallPackage(pkg, ver, gold_file="lib/libreadline.a")
    elif pkg == 'ncurses':
        success = InstallPackage(pkg, ver, gold_file="lib/libncurses.a")
    elif pkg == 'lua':
        readline_install = f"{install_dir}/{'readline'.upper()}/readline-" + \
                           f"{package_info['readline'][VERSION]}/chi_build"
        ncurses_install = f"{install_dir}/{'ncurses'.upper()}/ncurses-" + \
                          f"{package_info['ncurses'][VERSION]}/chi_build"
        success = InstallLuaPackage(pkg, ver, gold_file="lib/liblua.a",
                                    readline_install=readline_install,
                                    ncurses_install=ncurses_install)
    elif pkg == 'petsc':
        success = InstallPETSC(pkg, ver, gold_file="include/petsc")
    elif pkg == 'vtk':
        success = InstallVTK(pkg, ver, gold_file="include")
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
    print(error_beg +
          f"Failed to install {len(install_errors)} package(s)\n" +
          "Check the package's associated log file, i.e., PACKAGE/package_log.txt" +
          "for information as to why the install failed. You could also try to " +
          "manually install these packages")
    print(error_end)
    exit(1)

# ========================================== Create envvars file
module_file = open(f"{install_dir}/configure_deproots.sh", "w")

lua_version = package_info["lua"][VERSION]
petsc_version = package_info["petsc"][VERSION]
vtk_version = package_info["vtk"][VERSION]

lua_root = f"{install_dir}/LUA/lua-{lua_version}/install"
petsc_root = f"{install_dir}/PETSC/petsc-{petsc_version}-install"
vtk_dir = f"{install_dir}/VTK/vtk-{vtk_version}-install"

module_file.write(f'export LUA_ROOT={lua_root}\n')
module_file.write(f'export PETSC_ROOT={petsc_root}\n')
module_file.write(f'export VTK_DIR={vtk_dir}\n')

if os.path.exists(f"{vtk_dir}/lib64"):
    module_file.write(f'export LD_LIBRARY_PATH="{vtk_dir}/lib64":$LD_LIBRARY_PATH\n')
else:
    module_file.write(f'export LD_LIBRARY_PATH="{vtk_dir}/lib":$LD_LIBRARY_PATH\n')
module_file.write('echo "Enviroment set for compiling. ' +
                  'If recompiling changed sources execute"\n')
module_file.write('echo "     ./configure clean"\n')
module_file.write('echo " "\n')
module_file.write('echo "Otherwise just execute:"\n')
module_file.write('echo "     ./configure"\n')

module_file.close()

ExecSub("chmod u+x configure_deproots.sh", log_file)
log_file.close()

print("\n########## Chi-Tech Dependency install complete ##########")
print("\nWhen opening ChiTech in an IDE, the following environment variables" +
      " need to be set:\n")

print(f'LUA_ROOT = {lua_root}')
print(f'PETSC_ROOT = {petsc_root}')
print(f'VTK_DIR = {vtk_dir}')
print()
print(TextColors.WARNING +
      "When compiling ChiTech, in a terminal, the following environment variables" +
      " need to be set:\n" +
      TextColors.ENDC)

print(f'export LUA_ROOT="{lua_root}"')
print(f'export PETSC_ROOT="{petsc_root}"')
print(f'export VTK_DIR="{vtk_dir}"')
if os.path.exists(f"{vtk_dir}/lib64"):
    print(f'export LD_LIBRARY_PATH="{vtk_dir}/lib":$LD_LIBRARY_PATH')
else:
    print(f'export LD_LIBRARY_PATH="{vtk_dir}/lib64":$LD_LIBRARY_PATH')

print()
print(f"To set these terminal environment variables automatically, execute:")
print(f"    $. {install_dir}/configure_deproots.sh\n")
