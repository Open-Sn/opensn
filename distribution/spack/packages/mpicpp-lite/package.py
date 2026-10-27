from spack.package import *


class MpicppLite(CMakePackage):
    """Header-only MPI C++ wrappers used by OpenSn."""

    homepage = "https://github.com/andrsd/mpicpp-lite"
    git = "https://github.com/andrsd/mpicpp-lite.git"
    url = "https://github.com/andrsd/mpicpp-lite/archive/refs/tags/v3.1.0.tar.gz"

    version("3.1.0", tag="v3.1.0")

    depends_on("cmake@3.12:", type="build")
    depends_on("mpi")
