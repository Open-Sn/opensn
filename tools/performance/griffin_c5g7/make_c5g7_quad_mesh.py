#!/usr/bin/env python3
"""Write the 2D C5G7 all-quadrilateral mesh used by the Griffin GPU sweep paper.

N. Choi, J. Hanophy, and Y. Wang, "Development of a GPU-Based Discrete Ordinates
Transport Sweeper in the Griffin Reactor Physics Application," Nucl. Sci. Eng. (2026),
https://doi.org/10.1080/00295639.2026.2700698, Sec. IV.A and Fig. 4b.

Each 1.26 cm pin cell is meshed with the same topology as Fig. 4b:

* fuel pin (1156 assembly positions, including guide tubes and fission chambers):
  a 4x4 grid rotated 45 degrees (16 quads), a ring of 16 quads out to a 16-sided
  polygon approximating the r = 0.54 cm pin, and a ring of 16 moderator quads out
  to the pin-cell boundary (48 quads);
* moderator pin (1445 reflector positions): the rotated 4x4 grid and one ring of
  16 quads out to the pin-cell boundary (32 quads).

The total is 1156 * 48 + 1445 * 32 = 101,728 quads, matching the paper. The 16-gon
preserves the vertices of the r = 0.54 cm circle, so its area is 2.6% low; the paper
approximates the circle with flat sides in the same way.

Materials use the block IDs of test/python/.../transport_keigen/c5g7/c5g7.py:
0 water, 1 UO2, 2 7.0% MOX, 3 guide tube, 4 4.3% MOX, 5 8.7% MOX, 6 fission chamber.
The UO2 assembly occupies the corner at the origin, so xmin and ymin are reflecting.
"""

import argparse
import math

PITCH = 1.26
PIN_RADIUS = 0.54
NUM_ASSEMBLY_PINS = 17
NUM_CORE_PINS = 51

WATER, UO2, MOX70, GUIDE_TUBE, MOX43, MOX87, FISSION_CHAMBER = range(7)

GUIDE_TUBES = {(2, 5), (2, 8), (2, 11), (3, 3), (3, 13), (5, 2), (5, 5), (5, 8), (5, 11),
               (5, 14), (8, 2), (8, 5), (8, 11), (8, 14), (11, 2), (11, 5), (11, 8), (11, 11),
               (11, 14), (13, 3), (13, 13), (14, 5), (14, 8), (14, 11)}


def assembly_pin_material(kind, i, j):
    """Material of pin (i, j) in a 17x17 assembly of kind 'uo2' or 'mox'."""
    if (i, j) == (8, 8):
        return FISSION_CHAMBER
    if (i, j) in GUIDE_TUBES:
        return GUIDE_TUBE
    if kind == "uo2":
        return UO2
    if i in (0, 16) or j in (0, 16):
        return MOX43
    # The 8.7% zone is an octagon about the fission chamber.
    di, dj = abs(i - 8), abs(j - 8)
    return MOX87 if di <= 5 and dj <= 5 and di + dj <= 8 else MOX70


def pin_material(ix, iy):
    """Return the pin material, or None for a moderator (reflector) pin."""
    if ix >= 2 * NUM_ASSEMBLY_PINS or iy >= 2 * NUM_ASSEMBLY_PINS:
        return None
    kind = "uo2" if (ix < NUM_ASSEMBLY_PINS) == (iy < NUM_ASSEMBLY_PINS) else "mox"
    return assembly_pin_material(kind, ix % NUM_ASSEMBLY_PINS, iy % NUM_ASSEMBLY_PINS)


def ring_angles():
    return [2.0 * math.pi * k / 16.0 for k in range(16)]


def boundary_points(half):
    """16 points on the pin-cell boundary along the rays at k * 22.5 degrees."""
    pts = []
    for t in ring_angles():
        c, s = math.cos(t), math.sin(t)
        scale = half / max(abs(c), abs(s))
        pts.append((c * scale, s * scale))
    return pts


def diamond(a):
    """Rotated 4x4 grid with corners (+-a, 0), (0, +-a).

    Returns the 5x5 node grid and its 16 boundary nodes in counterclockwise order starting at
    (a, 0), so boundary node k lies near the ray at k * 22.5 degrees.
    """
    # Grid coordinates (u, v) in [0, 4]^2 map to the rotated square.
    def node(u, v):
        s, t = u / 4.0 - 0.5, v / 4.0 - 0.5
        return (a * (s - t), a * (s + t))

    grid = [[node(u, v) for v in range(5)] for u in range(5)]
    # Corner (u, v) = (4, 0) is (a, 0); walk the boundary counterclockwise.
    path = ([(4, v) for v in range(0, 4)] + [(u, 4) for u in range(4, 0, -1)]
            + [(0, v) for v in range(4, 0, -1)] + [(u, 0) for u in range(0, 4)])
    return grid, [grid[u][v] for u, v in path]


def pin_cells(is_fuel):
    """Return local (polygon, region) lists for one pin; region is 'inner' or 'outer'."""
    half = 0.5 * PITCH
    outer = boundary_points(half)
    cells = []
    if is_fuel:
        grid, ring0 = diamond(0.36)
        ring1 = [(PIN_RADIUS * math.cos(t), PIN_RADIUS * math.sin(t)) for t in ring_angles()]
        rings = [(ring0, ring1, "inner"), (ring1, outer, "outer")]
    else:
        grid, ring0 = diamond(0.42)
        rings = [(ring0, outer, "outer")]
    for u in range(4):
        for v in range(4):
            cells.append(([grid[u][v], grid[u + 1][v], grid[u + 1][v + 1], grid[u][v + 1]],
                          "inner" if is_fuel else "outer"))
    for inner_ring, outer_ring, region in rings:
        for k in range(16):
            k1 = (k + 1) % 16
            cells.append(([inner_ring[k], outer_ring[k], outer_ring[k1], inner_ring[k1]], region))
    return cells


def check_convex_ccw(poly):
    n = len(poly)
    for i in range(n):
        p0, p1, p2 = poly[i], poly[(i + 1) % n], poly[(i + 2) % n]
        cross = (p1[0] - p0[0]) * (p2[1] - p1[1]) - (p1[1] - p0[1]) * (p2[0] - p1[0])
        if cross <= 1.0e-12:
            raise RuntimeError(f"Non-convex or clockwise quad: {poly}")


def build_mesh():
    templates = {True: pin_cells(True), False: pin_cells(False)}
    for cells in templates.values():
        for poly, _ in cells:
            check_convex_ccw(poly)

    node_ids = {}
    nodes = []
    elements = []

    def node_id(x, y):
        key = (round(x, 9), round(y, 9))
        if key not in node_ids:
            node_ids[key] = len(nodes) + 1
            nodes.append((x, y))
        return node_ids[key]

    for iy in range(NUM_CORE_PINS):
        for ix in range(NUM_CORE_PINS):
            cx, cy = (ix + 0.5) * PITCH, (iy + 0.5) * PITCH
            mat = pin_material(ix, iy)
            for poly, region in templates[mat is not None]:
                cell_mat = mat if region == "inner" else WATER
                ids = [node_id(cx + px, cy + py) for px, py in poly]
                elements.append((cell_mat, ids))
    return nodes, elements


def write_gmsh22(path, nodes, elements):
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write(f"$Nodes\n{len(nodes)}\n")
        for i, (x, y) in enumerate(nodes, start=1):
            f.write(f"{i} {x:.12g} {y:.12g} 0\n")
        f.write("$EndNodes\n")
        f.write(f"$Elements\n{len(elements)}\n")
        for i, (mat, ids) in enumerate(elements, start=1):
            f.write(f"{i} 3 2 {mat} {mat} {' '.join(map(str, ids))}\n")
        f.write("$EndElements\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output", help="output Gmsh 2.2 .msh path")
    args = parser.parse_args()

    nodes, elements = build_mesh()
    write_gmsh22(args.output, nodes, elements)

    counts = [0] * 7
    for mat, _ in elements:
        counts[mat] += 1
    print(f"Wrote {args.output}: {len(nodes)} nodes, {len(elements)} quads")
    print("Quads per material (0..6):", counts)


if __name__ == "__main__":
    main()
