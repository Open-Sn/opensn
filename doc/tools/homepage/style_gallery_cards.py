# SPDX-License-Identifier: MIT
"""Build uniformly sized README gallery cards from committed source images."""

from pathlib import Path

import matplotlib
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps, ImageSequence


ROOT = Path(__file__).resolve().parents[3]
ASSETS = ROOT / "doc/assets"
BACKGROUND = "#000000"
DIVIDER = "#284052"
LINK = "#42d6c8"
CARD_SIZE = (1400, 1180)
PLOT_HEIGHT = 1010
MEDIA_ORIGIN = (100, 55)
MEDIA_SIZE = (1200, 900)
CONTENT_SIZE = (1120, 840)
BLACK = (0, 0, 0)

URBAN_RAMP = [
    "#104281",
    "#184f95",
    "#1c5cab",
    "#2a78d6",
    "#5598e7",
    "#86b6ef",
    "#cde2fb",
]


def font(size):
    """Use Matplotlib's bundled font for platform-independent rendering."""
    path = Path(matplotlib.get_data_path()) / "fonts/ttf/DejaVuSans-Bold.ttf"
    return ImageFont.truetype(str(path), size)


def exploration_band(image, label):
    """Replace only a card's caption area with one teal exploration label."""
    if image.size != CARD_SIZE:
        expected = f"{CARD_SIZE[0]}x{CARD_SIZE[1]}"
        raise ValueError(f"Expected a {expected} card, got {image.size}")
    image = image.convert("RGB")
    draw = ImageDraw.Draw(image)
    caption = (0, PLOT_HEIGHT, CARD_SIZE[0], CARD_SIZE[1])
    draw.rectangle(caption, fill=BACKGROUND)
    draw.line((74, PLOT_HEIGHT, 1326, PLOT_HEIGHT), fill=DIVIDER, width=2)
    draw.text((74, 1049), label, font=font(50), fill=LINK)
    return image


def trim_background(image, background=BLACK, tolerance=3):
    """Crop a uniform exterior background without changing interior pixels."""
    pixels = np.asarray(image.convert("RGB"), dtype=np.int16)
    color = np.asarray(background, dtype=np.int16)
    yy, xx = np.where(np.any(np.abs(pixels - color) > tolerance, axis=2))
    if len(xx) == 0:
        raise ValueError("Image contains no pixels distinct from its background")
    return image.crop((int(xx.min()), int(yy.min()), int(xx.max()) + 1, int(yy.max()) + 1))


def replace_exterior(image, source_color, target=BLACK):
    """Replace only the connected exterior background."""
    image = image.convert("RGB")
    corner = np.asarray(image.getpixel((0, 0)), dtype=np.int16)
    expected = np.asarray(source_color, dtype=np.int16)
    if np.any(np.abs(corner - expected) > 3):
        raise ValueError(
            f"Expected exterior color {source_color}, got {tuple(corner)}"
        )
    ImageDraw.floodfill(image, (0, 0), target, thresh=3)
    return image


def replace_background(image, source_color, target=BLACK, tolerance=3):
    """Replace a known solid plot background, including enclosed regions."""
    pixels = np.asarray(image.convert("RGB")).copy()
    source = np.asarray(source_color)
    mask = np.all(np.abs(pixels.astype(np.int16) - source) <= tolerance, axis=2)
    pixels[mask] = target
    return Image.fromarray(pixels, mode="RGB")


def style_transient_frame(image):
    """Blacken zero-value surfaces and remove the source wireframe."""
    pixels = np.asarray(image.convert("RGB")).copy()
    values = pixels.astype(np.int16)
    background = np.all(
        np.abs(values - np.asarray((24, 24, 28))) <= 3, axis=2
    )
    zero_surface = np.all(
        np.abs(values - np.asarray((52, 52, 60))) <= 3, axis=2
    )
    pixels[background | zero_surface] = BLACK

    spread = values.max(axis=2) - values.min(axis=2)
    brightness = values.mean(axis=2)
    wireframe = (spread <= 18) & (brightness >= 85)
    plot_region = np.zeros(wireframe.shape, dtype=bool)
    plot_region[75:640, 180:790] = True
    plot_region[:90, :425] = False
    wireframe &= plot_region
    median = np.asarray(
        Image.fromarray(pixels).filter(ImageFilter.MedianFilter(5))
    )
    pixels[wireframe] = median[wireframe]
    return Image.fromarray(pixels, mode="RGB")


def standard_card(image, label, *, trim=True):
    """Place one scientific image in the shared black media box."""
    image = image.convert("RGB")
    if trim:
        image = trim_background(image)
    image = ImageOps.contain(image, CONTENT_SIZE, Image.Resampling.LANCZOS)
    card = Image.new("RGB", CARD_SIZE, BACKGROUND)
    media = Image.new("RGB", MEDIA_SIZE, BLACK)
    media.paste(image, ((MEDIA_SIZE[0] - image.width) // 2,
                        (MEDIA_SIZE[1] - image.height) // 2))
    card.paste(media, MEDIA_ORIGIN)
    return exploration_band(card, label)


def c5g7_card():
    source = Image.open(ASSETS / "source-c5g7-panel.png")
    source = replace_exterior(source, (13, 27, 42))
    return standard_card(source, "Explore Reactor Physics", trim=False)


def wire_card():
    """Create a card from the OpenSn coiled-wire deposition solution."""
    source = Image.open(ASSETS / "source-wire-energy-deposition.png")
    source = replace_exterior(source, (13, 27, 42))
    return standard_card(source, "Explore Charged Particle Transport", trim=False)


def remap_urban_colors(image):
    """Map the committed example figure's documented blue ramp to inferno."""
    pixels = np.asarray(image.convert("RGB"), dtype=np.int32)
    flat = pixels.reshape(-1, 3)
    colors, inverse = np.unique(flat, axis=0, return_inverse=True)
    source_cmap = LinearSegmentedColormap.from_list("urban", URBAN_RAMP)
    positions = np.linspace(0.0, 1.0, 2048)
    source = np.rint(source_cmap(positions)[:, :3] * 255).astype(np.int32)
    target = np.rint(
        colormaps["inferno"](positions)[:, :3] * 255
    ).astype(np.uint8)
    nearest = np.empty(len(colors), dtype=np.int32)
    distance = np.empty(len(colors), dtype=np.int32)
    for begin in range(0, len(colors), 1024):
        end = min(begin + 1024, len(colors))
        difference = colors[begin:end, None, :] - source[None, :, :]
        squared = np.sum(difference * difference, axis=2)
        nearest[begin:end] = np.argmin(squared, axis=1)
        distance[begin:end] = np.min(squared, axis=1)
    mapped = colors.astype(np.uint8)
    ramp_pixels = distance <= 12**2
    mapped[ramp_pixels] = target[nearest[ramp_pixels]]
    return Image.fromarray(mapped[inverse].reshape(pixels.shape), mode="RGB")


def urban_card():
    """Create a card from the committed urban OpenSn dose-rate solution."""
    source = Image.open(ASSETS / "source-urban-dose-map.png").convert("RGB")
    opensn_panel = source.crop((0, 0, 840, source.height))
    opensn_panel = remap_urban_colors(opensn_panel)
    opensn_panel = replace_background(opensn_panel, (26, 26, 25))
    return standard_card(opensn_panel, "Explore Shielding and Streaming", trim=False)


def multiphysics_card():
    source = Image.open(ASSETS / "source-multiphysics-panel.png")
    source = replace_exterior(source, (13, 27, 42))
    return standard_card(source, "Explore Coupled Multiphysics", trim=False)


def transient_card():
    """Create a linked gallery animation with the source GIF's timing."""
    source = Image.open(ASSETS / "source-kobayashi-animation.gif")
    frames = []
    durations = []
    for source_frame in ImageSequence.Iterator(source):
        frame = style_transient_frame(source_frame)
        card = standard_card(frame,
                             "Explore Time Dependent Transport", trim=False)
        frames.append(card.convert("P", palette=Image.Palette.ADAPTIVE))
        durations.append(source_frame.info.get("duration", 120))
    output = ASSETS / "gallery-transient.gif"
    frames[0].save(
        output,
        save_all=True,
        append_images=frames[1:],
        duration=durations,
        loop=source.info.get("loop", 0),
        disposal=2,
        optimize=True,
    )


def tokamak_card():
    """Create a card from the documented OpenSn tokamak flux figure."""
    source = Image.open(ASSETS / "source-tokamak-flux.png").convert("RGB")
    source = replace_exterior(source, (13, 27, 41))
    composed = Image.new("RGB", MEDIA_SIZE, BLACK)
    title = ImageOps.contain(source.crop((250, 0, 1350, 100)),
                             (1000, 75), Image.Resampling.LANCZOS)
    geometry = ImageOps.contain(source.crop((150, 220, 1450, 1130)),
                                (1080, 650), Image.Resampling.LANCZOS)
    footer = ImageOps.contain(source.crop((40, 1140, 1120, 1400)),
                              (1040, 150), Image.Resampling.LANCZOS)
    composed.paste(title, ((MEDIA_SIZE[0] - title.width) // 2, 10))
    composed.paste(geometry, ((MEDIA_SIZE[0] - geometry.width) // 2, 85))
    composed.paste(footer, (60, 735))
    return standard_card(composed, "Explore Fusion Engineering", trim=False)


if __name__ == "__main__":
    c5g7_card().save(ASSETS / "gallery-c5g7.png", optimize=True)
    wire_card().save(ASSETS / "gallery-wire.png", optimize=True)
    urban_card().save(ASSETS / "gallery-urban.png", optimize=True)
    transient_card()
    multiphysics_card().save(ASSETS / "gallery-multiphysics.png", optimize=True)
    tokamak_card().save(ASSETS / "gallery-fusion.png", optimize=True)
