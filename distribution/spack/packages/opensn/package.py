from spack.package import *


class Opensn(CMakePackage):
    """OpenSn: a parallel scientific simulation engine for the linear Boltzmann equation."""

    homepage = "https://open-sn.github.io/opensn/"
    git = "https://github.com/Open-Sn/opensn.git"
    url = "https://github.com/Open-Sn/opensn/archive/refs/tags/v1.0.1.tar.gz"

    version("develop", branch="main")
    version("1.0.1", tag="v1.0.1")
    version("1.0.0", tag="v1.0.0")

    variant("python", default=True, description="Build the Python-enabled OpenSn console")
    variant("python_module", default=False, description="Build and install the pyopensn module")
    variant("cuda", default=False, description="Enable CUDA support")
    variant("hip", default=False, description="Enable HIP support")
    variant("sycl", default=False, description="Enable SYCL support")
    variant(
        "native",
        default=False,
        description="Build with OpenSn's native optimized CMake configuration",
    )
    variant(
        "sycl_compiler",
        default="icpx",
        values=("icpx", "acpp"),
        multi=False,
        description="SYCL compiler",
    )
    variant(
        "sycl_arch",
        default="spir64_gen",
        values=str,
        multi=False,
        description="SYCL target architecture passed to OpenSn CMake",
    )
    variant(
        "icpx_target_backend",
        default="none",
        values=str,
        multi=False,
        description="ICPX target backend option for AoT compilation, or 'none' to leave unset",
    )
    variant(
        "icpx_register_alloc_mode",
        default="none",
        values=str,
        multi=False,
        description="ICPX register allocation mode, or 'none' to leave unset",
    )
    variant(
        "acpp_platform",
        default="none",
        values=("none", "cuda", "rocm", "cpu"),
        multi=False,
        description="ACPP platform",
    )

    conflicts("+python_module", when="~python", msg="The pyopensn module requires +python")
    conflicts("+cuda", when="+hip", msg="CUDA, HIP, and SYCL support are mutually exclusive")
    conflicts("+cuda", when="+sycl", msg="CUDA, HIP, and SYCL support are mutually exclusive")
    conflicts("+hip", when="+sycl", msg="CUDA, HIP, and SYCL support are mutually exclusive")
    conflicts("+icpx_target_backend", when="+acpp_platform", msg="ICPX and ACPP are mutually exclusive")
    conflicts("+icpx_register_alloc_mode", when="+acpp_platform", msg="ICPX and ACPP are mutually exclusive")
    conflicts("+native", when="build_type=Debug", msg="+native sets CMAKE_BUILD_TYPE=Native")
    conflicts("+native", when="build_type=RelWithDebInfo", msg="+native sets CMAKE_BUILD_TYPE=Native")
    conflicts("+native", when="build_type=MinSizeRel", msg="+native sets CMAKE_BUILD_TYPE=Native")

    depends_on("cmake@3.29:", type="build")
    depends_on("gmake", type="build")
    depends_on("googletest@1.15:", type="build")
    depends_on("mpi")
    depends_on("mpicpp-lite@2.7:")
    depends_on("boost@1.86:")
    depends_on("hdf5@1.14:+hl+mpi")
    depends_on("vtk@9.3:+mpi io=exodusii")
    depends_on("caliper@2.13:+mpi~gotcha~kokkos")
    depends_on(
        "petsc@3.24.1:+mpi+double~complex+int64+hdf5+metis+ptscotch+hypre"
    )

    depends_on("cuda@12:", when="+cuda")
    depends_on("hip", when="+hip")
    depends_on("intel-oneapi-compilers", when="+sycl sycl_compiler=icpx", type="build")
    depends_on("adaptivecpp", when="+sycl sycl_compiler=acpp")

    depends_on("python@3.9:", when="+python")
    depends_on("py-pybind11", when="+python")
    depends_on("py-numpy", when="+python_module")
    depends_on("py-mpi4py", when="+python_module")

    test_source_paths = [
        "test/run_tests",
        "test/src",
        "test/python",
        "test/assets",
    ]

    def cmake_args(self):
        spec = self.spec
        args = [
            self.define_from_variant("OPENSN_WITH_PYTHON", "python"),
            self.define_from_variant("OPENSN_WITH_PYTHON_MODULE", "python_module"),
            self.define_from_variant("OPENSN_WITH_CUDA", "cuda"),
            self.define_from_variant("OPENSN_WITH_HIP", "hip"),
            self.define_from_variant("OPENSN_WITH_SYCL", "sycl"),
        ]

        if "+python" in spec:
            args.append(self.define("Python3_EXECUTABLE", spec["python"].command.path))

        if "+native" in spec:
            args.append(self.define("CMAKE_BUILD_TYPE", "Native"))

        if "+sycl" in spec:
            sycl_compiler = spec.variants["sycl_compiler"].value
            args.append(self.define("SYCL_COMPILER", sycl_compiler))
            args.append(self.define("SYCL_ARCHITECTURE", spec.variants["sycl_arch"].value))

            if sycl_compiler == "icpx":
                icpx_target_backend = spec.variants["icpx_target_backend"].value
                if icpx_target_backend != "none":
                    args.append(self.define("ICPX_TARGET_BACKEND", icpx_target_backend))
                icpx_register_alloc_mode = spec.variants["icpx_register_alloc_mode"].value
                if icpx_register_alloc_mode != "none":
                    args.append(self.define("ICPX_REGISTER_ALLOC_MODE", icpx_register_alloc_mode))
            elif sycl_compiler == "acpp":
                acpp_platform = spec.variants["acpp_platform"].value
                if acpp_platform != "none":
                    args.append(self.define("ACPP_PLATFORM", acpp_platform))

        return args

    def _unit_test_exe(self, installed=False):
        if installed:
            return join_path(self.prefix.libexec, "opensn", "opensn-unit")
        return join_path(self.build_directory, "test", "opensn-unit")

    def _opensn_exe(self, installed=False):
        if installed:
            return self.prefix.bin.opensn
        return join_path(self.build_directory, "python", "opensn")

    def _mpi_cmd(self):
        return "{0} -n ".format(join_path(self.spec["mpi"].prefix.bin, "mpiexec"))

    def _run_unit_tests(self, installed=False):
        unit_test = Executable(self._unit_test_exe(installed))
        with set_env(
            DYLD_LIBRARY_PATH=self.prefix.lib,
            LD_LIBRARY_PATH=self.prefix.lib,
        ):
            unit_test()

    def _run_regression_tests(self, test_root, installed=False):
        python = self.spec["python"].command
        python(
            join_path(test_root, "test", "run_tests"),
            "-d",
            join_path(test_root, "test", "python"),
            "--exe",
            self._opensn_exe(installed),
            "--mpi-cmd",
            self._mpi_cmd(),
            "-w3",
        )

    @run_after("build")
    @on_package_attributes(run_tests=True)
    def run_build_tests(self):
        if "~python" in self.spec:
            return
        self._run_unit_tests()
        self._run_regression_tests(self.stage.source_path)

    @run_after("install")
    def install_unit_test(self):
        if "~python" in self.spec:
            return
        mkdirp(join_path(self.prefix.libexec, "opensn"))
        install(self._unit_test_exe(), self._unit_test_exe(installed=True))

    @run_after("install")
    def cache_test_sources(self):
        if "~python" in self.spec:
            return
        cache_extra_test_sources(self, self.test_source_paths)

    def test_unit(self):
        """Run the OpenSn unit test executable."""
        if "~python" in self.spec:
            raise SkipTest("OpenSn unit tests are only built with +python")
        self._run_unit_tests(installed=True)

    def test_regression(self):
        """Run the OpenSn Python regression test suite."""
        if "~python" in self.spec:
            raise SkipTest("OpenSn regression tests require +python")
        self._run_regression_tests(self.test_suite.current_test_cache_dir, installed=True)
