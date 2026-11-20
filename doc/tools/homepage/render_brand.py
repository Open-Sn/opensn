# SPDX-License-Identifier: MIT
"""Generate the original OpenSn angular-sweep SVGs and transparent avatar.

The selected concept was explored with AI image generation, then reconstructed
as exact vector geometry. No simulation data or third-party logo is used.
"""

from pathlib import Path
import subprocess


DEST = Path(__file__).resolve().parents[3] / "doc/assets/brand"


def artwork(dark=False, mark=False):
    """Repeat one circular sector at three angles; share concentric mesh cuts."""
    navy = ("#f2f7fb", "#859fb7") if dark else ("#082a48", "#9bb1c6")
    teal = ("#42d6c8", "#238c90") if dark else ("#008c9f", "#96d5df")
    ink = "#f2f7fb" if dark else "#082a48"
    viewbox = "0 0 200 200" if mark else "0 0 620 112"
    transform = "translate(32 165) scale(.95)" if mark else "translate(14 91) scale(.58)"
    wordmark = "" if mark else (
        '<text x="127" y="82" font-family="DejaVu Sans, Arial, sans-serif" '
        f'font-size="105" font-weight="700" letter-spacing="-3" fill="{ink}">OpenSn</text>'
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg"
  xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="{viewbox}"
  role="img" aria-labelledby="title desc">
  <title id="title">OpenSn angular sweep</title>
  <desc id="desc">Three identical angular sectors extend from an orange particle,
  with three concentric mesh divisions and fading intensity.</desc>
  <!-- SPDX-License-Identifier: MIT. Original artwork designed for OpenSn. -->
  <!-- Radius 160; half-angle atan(1/8); axes -50, -30, -10 degrees.
       Mesh cuts at radii 40, 80, 120. All sectors share one geometry. -->
  <defs>
    <linearGradient id="navy" gradientUnits="userSpaceOnUse" x1="0" y1="0" x2="160" y2="0">
      <stop offset="0" stop-color="{navy[0]}"/>
      <stop offset="1" stop-color="{navy[1]}"/>
    </linearGradient>
    <linearGradient id="teal" gradientUnits="userSpaceOnUse" x1="0" y1="0" x2="160" y2="0">
      <stop offset="0" stop-color="{teal[0]}"/>
      <stop offset="1" stop-color="{teal[1]}"/>
    </linearGradient>
    <mask id="mesh" maskUnits="userSpaceOnUse" x="-2" y="-25"
      width="170" height="50" style="mask-type:luminance">
      <rect x="-2" y="-25" width="170" height="50" fill="white"/>
      <g fill="none" stroke="black" stroke-width="3">
        <circle r="40"/><circle r="80"/><circle r="120"/>
      </g>
    </mask>
    <g id="cone" mask="url(#mesh)">
      <path d="M0 0 158.764460274 -19.845557534 A160 160 0 0 1 158.764460274 19.845557534Z"/>
    </g>
  </defs>
  <g transform="{transform}">
    <use xlink:href="#cone" transform="rotate(-50)" fill="url(#navy)"/>
    <use xlink:href="#cone" transform="rotate(-30)" fill="url(#teal)"/>
    <use xlink:href="#cone" transform="rotate(-10)" fill="url(#navy)"/>
    <circle r="12" fill="#ff8800"/>
  </g>
  {wordmark}
</svg>
'''


if __name__ == "__main__":
    for filename, dark, mark in (
        ("logo-light.svg", False, False),
        ("logo-dark.svg", True, False),
        ("mark.svg", False, True),
    ):
        (DEST / filename).write_text(artwork(dark=dark, mark=mark), encoding="utf-8")
    subprocess.run([
        "rsvg-convert", "-w", "512", "-h", "512",
        "-o", str(DEST / "mark.png"), str(DEST / "mark.svg"),
    ], check=True)
