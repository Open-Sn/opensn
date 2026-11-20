#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A split-mesh file that cannot be written raises an error (1 rank, 2 parts).

The process file-size limit is set below the size of the first part file, so writing that file
fails with EFBIG. File streams do not throw, so SplitFileMeshGenerator must check the stream after
it closes the file. Without that check the generator reports success, and with fewer ranks than
parts the run ends with a successful exit status and an incomplete split mesh.
"""

import os
import resource
import shutil
import signal
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator, SplitFileMeshGenerator

if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Incorrect number of processors. Expected 1 processor but got {size}.")

    # 20^3 hexahedra give part files of about 1.5 MiB each, above the 1 MiB limit. The log is much
    # smaller than the limit.
    num_cells = 20
    nodes = [i / num_cells for i in range(num_cells + 1)]
    file_size_limit = 1024 * 1024
    # The directory is relative to the working directory, so the test writes nothing elsewhere.
    split_mesh_dir = "split_mesh_write_error"

    # Exceeding the limit raises SIGXFSZ, which terminates the process unless it is ignored. The
    # write then fails with EFBIG instead.
    old_handler = signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
    old_limit = resource.getrlimit(resource.RLIMIT_FSIZE)
    resource.setrlimit(resource.RLIMIT_FSIZE, (file_size_limit, old_limit[1]))
    try:
        SplitFileMeshGenerator(
            inputs=[OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])],
            num_partitions=2,
            split_mesh_dir_path=split_mesh_dir,
        ).Execute()
        raised = 0
    except RuntimeError as error:
        raised = 1 if "Failed to write" in str(error) else 0
    finally:
        resource.setrlimit(resource.RLIMIT_FSIZE, old_limit)
        signal.signal(signal.SIGXFSZ, old_handler)
        # Remove the incomplete part files.
        shutil.rmtree(split_mesh_dir, ignore_errors=True)

    print(f"Split mesh write failure raised = {raised}")
