# Homepage images

This directory contains the images used by the repository `README.md`. The
small rendering scripts are in `doc/tools/homepage`.

Run the commands below from the repository root.

## Requirements

The scripts require Python with NumPy, Pillow, Matplotlib, PyVista, and VTK.
Generating the SVG logo and hero image also requires `rsvg-convert` from
librsvg. Recomputing the C5G7 image requires a built OpenSn console.

## Source images

Four panels are made from figures committed to OpenSnExamples at revision
`7e2d56a39`:

| Local source | OpenSnExamples source |
| --- | --- |
| `source-wire-energy-deposition.png` | `Coiled_Wire_Electron_Photon/assets/energy-deposition-dark.png` |
| `source-urban-dose-map.png` | `Urban_Source/assets/dose-map-dark.png` |
| `source-kobayashi-animation.gif` | `Kobayashi_Time_Dependent/kobayashi_td_animation.gif` |
| `source-tokamak-flux.png` | `Tokamak_Fusion_Engineering/assets/quarter-flux-dark.png` |

The hero image also uses the tokamak figure. The OpenSnExamples material is
distributed under the MIT license reproduced in `OpenSnExamples-LICENSE`.

`source-c5g7-panel.png` and `source-multiphysics-panel.png` are OpenSn field
renderings written by `render_gallery.py` and `render_fields.py`. They are
committed so that the cards can be rebuilt without rerunning OpenSn.

## Generate the logo

```sh
python doc/tools/homepage/render_brand.py
```

This writes the light and dark wordmarks, SVG mark, and PNG mark under
`doc/assets/brand`.

## Generate the example panels

```sh
python doc/tools/homepage/style_gallery_cards.py
```

This builds all six cards from the committed source images.

All six cards are 1400 by 1180 pixels and place their scientific figure in the
same 1200 by 900 pixel media area. The original scales, units, legends, and
annotations are retained while the figures are fitted to comparable bounds.

The script preserves the Kobayashi animation timing and loop behavior. It
removes the wireframe box from each frame and draws the gray surfaces, which
show values below the color-scale minimum, in black. For the
urban panel, it selects the OpenSn half of the comparison figure and converts
the documented blue color ramp to inferno without changing the normalized
color values, axes, geometry, source marker, scale, or units.

## Regenerate C5G7 and the hero image

Export the C5G7 thermal-group scalar flux:

```sh
build/python/opensn -i doc/tools/homepage/export_gallery.py
```

Then render the C5G7 source panel and the hero, and rebuild the cards. The
hero places the C5G7 field and the
tokamak geometry from `source-tokamak-flux.png` on one shared black two-panel
canvas:

```sh
python doc/tools/homepage/render_gallery.py \
  --data-directory gallery-output
python doc/tools/homepage/style_gallery_cards.py
```

The C5G7 field is normalized by its global maximum for display. The renderer
rejects non-finite or negative values.

## Regenerate the multiphysics panel

Given the energy-deposition PVTU file from the OpenFOAM electron-converter
example:

```sh
python doc/tools/homepage/render_fields.py \
  --converter /path/to/showcase_energy_deposition.pvtu
python doc/tools/homepage/style_gallery_cards.py
```

The image shows linear-tetrahedron cell-average energy deposition normalized
by its maximum. Values below the logarithmic display floor use the lowest
inferno color. The underlying OpenSn output is not modified or smoothed.
