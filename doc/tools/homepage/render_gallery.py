# SPDX-License-Identifier: MIT
"""Render README artwork from actual OpenSn outputs, never synthetic fields.

Requires numpy, Pillow, matplotlib, PyVista/VTK, and rsvg-convert. Raw meshes and
solution dumps are deliberately not committed. See doc/assets/README.md.
"""

import argparse
from pathlib import Path
import subprocess

from matplotlib import colormaps
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parents[3]
DEST = ROOT / "doc/assets"
BACKGROUND = "#0d1b2a"
INK = "#f2f7fb"
MUTED = "#a4b7c8"


def font(size, bold=False):
    """Use matplotlib's bundled fonts so rendering does not depend on system fonts."""
    import matplotlib
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    return ImageFont.truetype(str(Path(matplotlib.get_data_path()) / "fonts/ttf" / name), size)


def plotter():
    import pyvista as pv

    scene = pv.Plotter(off_screen=True, window_size=(1400, 1000))
    scene.set_background(BACKGROUND)
    scene.enable_anti_aliasing("ssaa")
    return scene


def screenshot(scene):
    pixels = scene.screenshot(return_img=True)
    scene.close()
    return Image.fromarray(pixels).convert("RGB")


def color_scale(image, cmap="inferno", label="Relative cell-average field"):
    drawing = ImageDraw.Draw(image)
    x, y, width = 74, 917, 510
    colors = colormaps[cmap](np.linspace(0, 1, width))[:, :3]
    strip = Image.fromarray((colors[None, :, :] * 255).astype(np.uint8)).resize((width, 12))
    image.paste(strip, (x, y))
    drawing.text((x, y - 30), label, font=font(19), fill=MUTED)
    drawing.text((x, y + 18), "0", font=font(18), fill=MUTED)
    drawing.text((x + width - 10, y + 18), "1", font=font(18), fill=MUTED)


def trim_plot(image, padding=12):
    """Remove empty plot margins, preserving every visible data pixel."""
    pixels = np.asarray(image).astype(int)
    background = np.array([13, 27, 42])
    yy, xx = np.where(np.any(np.abs(pixels - background) > 3, axis=2))
    return image.crop((max(0, int(xx.min()) - padding), max(0, int(yy.min()) - padding),
                       min(image.width, int(xx.max()) + padding + 1),
                       min(image.height, int(yy.max()) + padding + 1)))


def replace_background(image, source_color, tolerance=3):
    """Replace a known solid plot background, including enclosed regions."""
    pixels = np.asarray(image.convert("RGB")).copy()
    source = np.asarray(source_color)
    mask = np.all(np.abs(pixels.astype(np.int16) - source) <= tolerance, axis=2)
    pixels[mask] = (0, 0, 0)
    return Image.fromarray(pixels, mode="RGB")


def compose_hero(c5g7_image):
    """Build the hero with a shared C5G7/tokamak two-panel canvas."""
    hero = Image.new("RGB", (1800, 680), BACKGROUND)
    drawing = ImageDraw.Draw(hero)

    # Ease the navy text field into the black visualization field without
    # tinting the scientific images themselves.
    gradient_start = 550
    visualization_start = 850
    navy = np.asarray((13, 27, 42), dtype=float)
    for x in range(gradient_start, visualization_start):
        fraction = (x - gradient_start) / (visualization_start - gradient_start - 1)
        color = tuple(np.rint(navy * (1.0 - fraction)).astype(np.uint8))
        drawing.line((x, 0, x, hero.height), fill=color)
    drawing.rectangle((visualization_start, 0, hero.width, hero.height), fill="black")

    logo = DEST / "brand/logo-dark.svg"
    raster = DEST / ".logo-dark-homepage.png"
    subprocess.run(["rsvg-convert", "-w", "620", "-o", str(raster), str(logo)], check=True)
    wordmark = Image.open(raster).convert("RGBA")
    hero.paste(wordmark, (75, 60), wordmark)
    raster.unlink()
    drawing.text((84, 225), "Radiation transport.", font=font(64, True), fill=INK)
    drawing.text((84, 309), "Across every scale.", font=font(64, True), fill=INK)
    drawing.text((88, 425), "Deterministic physics. Parallel by design.", font=font(25), fill=MUTED)
    drawing.line((88, 506, 660, 506), fill="#284052", width=2)
    drawing.text((88, 537), "C++  /  Python  /  CPUs & GPUs", font=font(22), fill="#42d6c8")

    def floating_panel(image, y_offset=0):
        image = ImageOps.contain(image, (360, 360), Image.Resampling.LANCZOS)
        panel = Image.new("RGB", (400, 400), "black")
        panel.paste(image, ((400 - image.width) // 2,
                            (400 - image.height) // 2 + y_offset))
        return panel

    fields = Image.new("RGB", (860, 400), "black")
    c5g7_panel = floating_panel(c5g7_image, y_offset=-15)
    fields.paste(c5g7_panel, (0, 0))

    # Keep only the tokamak geometry. The title, annotations, and color bar of
    # the source figure are unreadable at panel size. Replace the source canvas,
    # including the region enclosed by the coils, so both panels share black.
    tokamak = Image.open(DEST / "source-tokamak-flux.png").convert("RGB")
    tokamak = tokamak.crop((0, 100, tokamak.width, 1150))
    tokamak = trim_plot(tokamak, padding=0)
    tokamak = replace_background(tokamak, (13, 27, 41))
    tokamak_panel = floating_panel(tokamak, y_offset=16)
    fields.paste(tokamak_panel, (460, 0))

    fields_x = 850 + (hero.width - 850 - fields.width) // 2
    hero.paste(fields, (fields_x, 170))
    drawing.text((fields_x + 200, 602), "C5G7 Thermal Flux",
                 font=font(20), fill=MUTED,
                 anchor="mt")
    drawing.text((fields_x + 660, 602), "Tokamak Neutron Flux",
                 font=font(20), fill=MUTED,
                 anchor="mt")
    hero.save(DEST / "hero.png", optimize=True)


def render(args):
    import pyvista as pv

    DEST.mkdir(parents=True, exist_ok=True)

    core_file = args.data_directory / "c5g7_thermal.pvtu"
    core = pv.read(core_file)
    fields = [name for name in core.point_data
              if name not in {"Block", "Partition", "vtkGhostType"}]
    if len(fields) != 1:
        raise ValueError(f"Expected only the exported thermal-group field, got {fields}.")
    field = fields[0]
    values = core.point_data[field]
    if not np.all(np.isfinite(values)) or np.any(values < 0) or values.max() <= 0:
        raise ValueError("C5G7 contains invalid or negative flux; do not hide it in the image.")
    core.point_data["relative"] = values / values.max()
    scene = plotter()
    scene.add_mesh(core, scalars="relative", preference="point", cmap="inferno",
                   clim=(0, 1), show_scalar_bar=False, lighting=False)
    scene.view_xy()
    scene.enable_parallel_projection()
    scene.camera.zoom(1.08)
    core_image = screenshot(scene)
    core_raw = core_image.copy()
    color_scale(core_image, label="Relative scalar flux · thermal group")
    core_image.save(DEST / "source-c5g7-panel.png", optimize=True)

    compose_hero(trim_plot(core_raw, padding=0))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-directory", type=Path, required=True)
    render(parser.parse_args())
