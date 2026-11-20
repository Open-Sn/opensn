# SPDX-License-Identifier: MIT
"""Render actual OpenSn gallery fields with a shared inferno palette.

Run with a Python environment containing PyVista/VTK, NumPy, Pillow, Matplotlib.
No synthetic fields, smoothing, or absolute values. Floors affect display only.
"""
import argparse
from pathlib import Path

import matplotlib
from matplotlib import colormaps
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pyvista as pv


ROOT = Path(__file__).resolve().parents[3]
ASSETS = ROOT / 'doc/assets'
BACKGROUND = '#0d1b2a'
INK = '#f2f7fb'
MUTED = '#a4b7c8'


def font(size, bold=False):
    name = 'DejaVuSans-Bold.ttf' if bold else 'DejaVuSans.ttf'
    return ImageFont.truetype(str(Path(matplotlib.get_data_path()) / 'fonts/ttf' / name), size)


def prepare(path, field, floor):
    mesh = pv.read(path)
    if not isinstance(mesh, pv.UnstructuredGrid):
        raise ValueError('Expected one OpenSn unstructured grid')
    if field not in mesh.cell_data:
        if set(mesh.celltypes) != {pv.CellType.TETRA}:
            raise ValueError('Exact nodal averaging is justified only for linear tetrahedra')
        mesh = mesh.point_data_to_cell_data(pass_point_data=False)
    values = np.asarray(mesh.cell_data[field])
    if not np.all(np.isfinite(values)) or values.max() <= 0:
        raise ValueError('Nonfinite field or no positive values')
    maximum = float(values.max())
    relative = values / maximum
    display = np.maximum(relative, floor)  # Visualization only; raw field is retained.
    mesh.cell_data['relative'] = display
    # OpenSn's discontinuous export repeats cell vertices. Merge coincident
    # geometry points only so interior tetrahedral faces are not drawn as surfaces.
    count = mesh.n_cells
    mesh = mesh.clean(tolerance=0)
    if mesh.n_cells != count:
        raise ValueError('Geometry cleaning changed the cell count')
    if not np.array_equal(mesh.cell_data[field], values):
        raise ValueError('Geometry cleaning changed scientific cell values')
    return mesh, (floor, 1.0)


def scene():
    plot = pv.Plotter(off_screen=True, window_size=(1400, 1000))
    plot.set_background(BACKGROUND)
    plot.enable_anti_aliasing('ssaa')
    return plot


def add_field(plot, mesh, limits, **kwargs):
    plot.add_mesh(mesh, scalars='relative', preference='cell', cmap='inferno',
                  clim=limits, log_scale=True, nan_color='#59616b',
                  show_scalar_bar=False, lighting=False, **kwargs)


def save_source(plot, limits, output, left_image=None):
    image = Image.fromarray(plot.screenshot(return_img=True)).convert('RGB')
    plot.close()
    if left_image is not None:
        inset = left_image.resize((350, 250), Image.Resampling.LANCZOS)
        image.paste(inset, (60, 105))
    draw = ImageDraw.Draw(image)
    if left_image is not None:
        draw.rectangle((59, 69, 411, 356), outline='#284052', width=2)
        draw.text((74, 80), 'Geometry', font=font(19, True), fill=INK)
        draw.text((500, 55), 'Radiation energy deposition',
                  font=font(26, True), fill=INK)
    strip = colormaps['inferno'](np.linspace(0, 1, 510))[:, :3]
    strip = Image.fromarray((strip[None, :, :] * 255).astype(np.uint8)).resize((510, 12))
    image.paste(strip, (74, 933))
    draw.text((74, 900), 'Relative cell-average field · log scale',
              font=font(19), fill=MUTED)
    draw.text((74, 952), f'{limits[0]:.1e}', font=font(18), fill=MUTED)
    draw.text((574, 952), '1', font=font(18), fill=MUTED)
    draw.text((750, 933), f'Display floor: {limits[0]:.0e}',
              font=font(18), fill=MUTED)
    image.save(output, optimize=True)


def converter(path):
    mesh, limits = prepare(path, 'waterblock_energy_deposition', floor=1e-8)
    plot = scene()
    for block, color, opacity in ((101, '#f0b46b', 1), (102, '#a4b7c8', 0.2),
                                  (106, '#42d6c8', 1), (105, '#6c8295', 0.2)):
        part = mesh.threshold([block, block], scalars='Block')
        plot.add_mesh(part, color=color, opacity=opacity, ambient=0.45,
                      diffuse=0.55, show_edges=False)
    plot.camera_position = [(-0.065, -0.075, 0.065), (0.019, 0, 0), (0, 0, 1)]
    plot.enable_parallel_projection()
    plot.reset_camera()
    plot.camera.zoom(0.9)
    geometry = Image.fromarray(plot.screenshot(return_img=True)).convert('RGB')
    plot.close()
    plot = scene()
    # Expose the target and converter interior without altering their computed values.
    for block in (101, 102):
        part = mesh.threshold([block, block], scalars='Block')
        add_field(plot, part.clip(normal='z', origin=(0, 0, 0), invert=True), limits)
    water = mesh.threshold([106, 106], scalars='Block')
    plot.add_mesh(water, color='#42d6c8', opacity=0.8, show_edges=False,
                  ambient=0.5, diffuse=0.5)
    housing = mesh.threshold([103, 103], scalars='Block')
    plot.add_mesh(housing.extract_surface().extract_feature_edges(), color='#6c8295',
                  opacity=0.35, line_width=1)
    shield = mesh.threshold([105, 105], scalars='Block')
    plot.add_mesh(shield, color='#8798a8', opacity=0.15, show_edges=False)
    plot.camera_position = [(-0.065, -0.075, 0.065), (0.019, 0, 0), (0, 0, 1)]
    plot.enable_parallel_projection()
    plot.reset_camera()
    plot.camera.zoom(0.95)
    save_source(plot, limits, ASSETS / 'source-multiphysics-panel.png', left_image=geometry)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--converter', type=Path, required=True)
    converter(parser.parse_args().converter)
