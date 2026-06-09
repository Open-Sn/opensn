#!/usr/bin/env python3

from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[5]
MESH_FILE = REPO_ROOT / "test/assets/mesh/uncollided_cube_coarse.msh"
UNCOLLIDED_FILE = SCRIPT_DIR / "uncollided.h5"

SOURCE_LOCATION = (0.010, 0.012, 0.014)
SAMPLE_POINT = (0.024, 0.016, 0.008)
DOMAIN_BOUNDS = (-0.001, 0.033)


def remove_file(path):
    try:
        Path(path).unlink()
    except FileNotFoundError:
        pass


def make_mesh(FromFileMeshGenerator):
    grid = FromFileMeshGenerator(filename=str(MESH_FILE)).Execute()
    grid.SetUniformBlockID(0)
    return grid


def make_xs(MultiGroupXS):
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=0.8, c=0.45)
    return xs


def make_whole_domain(RPPLogicalVolume):
    xmin, xmax = DOMAIN_BOUNDS
    return RPPLogicalVolume(
        xmin=xmin,
        xmax=xmax,
        ymin=xmin,
        ymax=xmax,
        zmin=xmin,
        zmax=xmax,
    )
