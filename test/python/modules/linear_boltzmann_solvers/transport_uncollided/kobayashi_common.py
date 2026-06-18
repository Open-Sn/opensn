"""Common OpenSn construction helpers for Kobayashi Problem 3."""

import os
from pathlib import Path

import kobayashi_dog_leg_model as model

MESH_DIR = Path(__file__).resolve().parents[4] / "assets" / "mesh" / "kobayashi"


def mesh_path(mesh_filename=None):
    file_name = mesh_filename or os.environ.get("KOBAYASHI_MESH_FILE", model.MESH_FILENAME)
    return MESH_DIR / file_name


def make_mesh(FromFileMeshGenerator):
    path = mesh_path()
    if not path.exists():
        raise RuntimeError(f"Missing {path}. Generate it from kobayashi_dog_leg.geo.")
    grid = FromFileMeshGenerator(filename=str(path)).Execute()
    grid.SetOrthogonalBoundaries()
    return grid


def make_cross_sections(case, MultiGroupXS):
    result = {}
    for name, (sigma_t, scattering_ratio) in model.cross_sections(case).items():
        xs = MultiGroupXS()
        xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=scattering_ratio)
        result[name] = xs
    return result


def make_xs_map(cross_sections):
    return [
        {"block_ids": [block_id], "xs": cross_sections[name]}
        for name, block_id in model.MATERIAL_IDS.items()
    ]
