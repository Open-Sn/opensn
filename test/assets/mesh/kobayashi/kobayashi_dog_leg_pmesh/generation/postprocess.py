#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""Clip a padded pMesh Voronoi mesh to the exact Kobayashi Problem 3 regions."""

import argparse
import math
import re
from collections import defaultdict
from pathlib import Path

TOL = 1.0e-6
# Material-interface planes and their active transverse extents. Applying a
# plane to any cell overlapping the corresponding rectangle preserves exact
# interfaces without slicing unrelated shield cells.
MATERIAL_PLANES = ((0, 10.0), (0, 30.0), (0, 40.0),
                   (1, 10.0), (1, 50.0), (1, 60.0),
                   (2, 10.0), (2, 30.0), (2, 40.0))
OUTER_PLANES = ((0, 0.0), (0, 60.0), (1, 0.0),
                (1, 100.0), (2, 0.0), (2, 60.0))
CUTS = tuple((axis, value, None) for axis, value in OUTER_PLANES + MATERIAL_PLANES)


def read_list(path):
    text = path.read_text()
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.S)
    text = re.sub(r'//.*', '', text)
    matches = list(re.finditer(r'(?m)^\s*(\d+)\s*\n\s*\(', text))
    if not matches:
        raise RuntimeError(f'Cannot find OpenFOAM list in {path}')
    match = matches[-1]
    return int(match.group(1)), text[match.end():text.rfind(')')]


def load_mesh(path):
    _, point_text = read_list(path / 'points')
    points = [tuple(map(float, item.split()))
              for item in re.findall(r'\(([^()]*)\)', point_text)
              if len(item.split()) == 3]
    _, face_text = read_list(path / 'faces')
    faces = [list(map(int, item.split()))
             for item in re.findall(r'\d+\(([^()]*)\)', face_text)]
    _, owner_text = read_list(path / 'owner')
    owners = list(map(int, re.findall(r'-?\d+', owner_text)))
    _, neighbour_text = read_list(path / 'neighbour')
    neighbours = [value for value in map(int, re.findall(r'-?\d+', neighbour_text))
                  if value >= 0]
    points = [list(point) for point in points]
    plane_by_name = {'xmin': (0, 0.0), 'xmax': (0, 60.0),
                     'ymin': (1, 0.0), 'ymax': (1, 100.0),
                     'zmin': (2, 0.0), 'zmax': (2, 60.0)}
    boundary_text = re.sub(r'/\*.*?\*/', '', (path / 'boundary').read_text(), flags=re.S)
    for name, block in re.findall(r'(\w+)\s*\{(.*?)\}', boundary_text, flags=re.S):
        if name not in plane_by_name:
            continue
        count = int(re.search(r'nFaces\s+(\d+)', block).group(1))
        start = int(re.search(r'startFace\s+(\d+)', block).group(1))
        axis, value = plane_by_name[name]
        for face in faces[start:start + count]:
            for vertex in face:
                points[vertex][axis] = value
    points = [tuple(point) for point in points]
    if len(owners) != len(faces):
        raise RuntimeError('owner/face count mismatch')
    ncells = max(owners + neighbours) + 1
    cells = [[] for _ in range(ncells)]
    for face_id, face in enumerate(faces):
        cells[owners[face_id]].append([points[v] for v in face])
        if face_id < len(neighbours):
            cells[neighbours[face_id]].append([points[v] for v in reversed(face)])
    return cells


def close(a, b):
    return all(abs(a[i] - b[i]) <= TOL for i in range(3))


def clean_polygon(poly):
    result = []
    for point in poly:
        point = tuple(0.0 if abs(v) <= TOL else v for v in point)
        if not result or not close(point, result[-1]):
            result.append(point)
    if len(result) > 1 and close(result[0], result[-1]):
        result.pop()
    changed = True
    while changed and len(result) > 3:
        changed = False
        for i in range(len(result)):
            a, b, c = result[i - 1], result[i], result[(i + 1) % len(result)]
            ab = tuple(b[j] - a[j] for j in range(3))
            bc = tuple(c[j] - b[j] for j in range(3))
            cross = (ab[1] * bc[2] - ab[2] * bc[1],
                     ab[2] * bc[0] - ab[0] * bc[2],
                     ab[0] * bc[1] - ab[1] * bc[0])
            cross_norm = math.sqrt(sum(v * v for v in cross))
            scale = math.sqrt(sum(v * v for v in ab) * sum(v * v for v in bc))
            if cross_norm <= TOL * max(scale, 1.0):
                result.pop(i)
                changed = True
                break
    return result if len(result) >= 3 else []


def clipped_face(face, axis, value, keep_low):
    result = []
    intersections = []
    for start, end in zip(face, face[1:] + face[:1]):
        ds = start[axis] - value
        de = end[axis] - value
        if abs(ds) <= TOL:
            ds = 0.0
        if abs(de) <= TOL:
            de = 0.0
        start_in = ds <= 0.0 if keep_low else ds >= 0.0
        end_in = de <= 0.0 if keep_low else de >= 0.0
        if start_in:
            result.append(start)
        if start_in != end_in:
            fraction = (value - start[axis]) / (end[axis] - start[axis])
            point = tuple(value if i == axis else
                          start[i] + fraction * (end[i] - start[i])
                          for i in range(3))
            result.append(point)
            intersections.append(point)
    return clean_polygon(result), intersections


def cap_polygon(points, axis, outward_positive):
    unique = {}
    for point in points:
        key = tuple(round(v / TOL) for v in point)
        unique[key] = point
    cap = list(unique.values())
    if len(cap) < 3:
        return []
    basis = ((1, 2), (2, 0), (0, 1))[axis]
    centre = tuple(sum(p[i] for p in cap) / len(cap) for i in range(3))
    cap.sort(
        key=lambda p: math.atan2(
            p[basis[1]] - centre[basis[1]], p[basis[0]] - centre[basis[0]]
        )
    )
    if not outward_positive:
        cap.reverse()
    return clean_polygon(cap)


def clip_cell(cell, axis, value, keep_low):
    faces = []
    intersections = []
    for face in cell:
        clipped, face_intersections = clipped_face(face, axis, value, keep_low)
        if clipped:
            faces.append(clipped)
        intersections.extend(face_intersections)
    cap = cap_polygon(intersections, axis, outward_positive=keep_low)
    if cap:
        faces.append(cap)
    return faces


def split_cell(cell, axis, value, transverse_ranges):
    if transverse_ranges is not None:
        transverse_axes = [candidate for candidate in range(3) if candidate != axis]
        for other_axis, bounds in zip(transverse_axes, transverse_ranges):
            coordinates = [point[other_axis] for face in cell for point in face]
            if max(coordinates) <= bounds[0] + TOL or min(coordinates) >= bounds[1] - TOL:
                return [cell]
    coordinates = [point[axis] for face in cell for point in face]
    if max(coordinates) <= value + TOL or min(coordinates) >= value - TOL:
        return [cell]
    low = clip_cell(cell, axis, value, True)
    high = clip_cell(cell, axis, value, False)
    return [part for part in (low, high) if part]


def signed_volume(cell):
    volume = 0.0
    for face in cell:
        a = face[0]
        for index in range(1, len(face) - 1):
            b, c = face[index], face[index + 1]
            cross = (b[1] * c[2] - b[2] * c[1],
                     b[2] * c[0] - b[0] * c[2],
                     b[0] * c[1] - b[1] * c[0])
            volume += a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]
    return volume / 6.0


def material_id(cell):
    unique = {}
    for face in cell:
        for point in face:
            unique[tuple(round(v / TOL) for v in point)] = point
    centre = tuple(sum(p[i] for p in unique.values()) / len(unique) for i in range(3))
    x, y, z = centre
    if x < 10.0 - TOL and y < 10.0 - TOL and z < 10.0 - TOL:
        return 1
    in_duct = ((x < 10.0 - TOL and y < 50.0 - TOL and z < 10.0 - TOL)
               or (x < 30.0 - TOL and 50.0 + TOL < y < 60.0 - TOL and z < 10.0 - TOL)
               or (30.0 + TOL < x < 40.0 - TOL and 50.0 + TOL < y < 60.0 - TOL
                   and z < 30.0 - TOL)
               or (30.0 + TOL < x < 40.0 - TOL and 50.0 + TOL < y
                   and 30.0 + TOL < z < 40.0 - TOL))
    return 2 if in_duct else 3


def polygon_area_2d(poly):
    return 0.5 * sum(poly[i][0] * poly[(i + 1) % len(poly)][1]
                     - poly[(i + 1) % len(poly)][0] * poly[i][1]
                     for i in range(len(poly)))


def intersect_convex_2d(subject, clip):
    if polygon_area_2d(clip) < 0.0:
        clip = list(reversed(clip))
    output = subject
    for a, b in zip(clip, clip[1:] + clip[:1]):
        input_poly, output = output, []
        if not input_poly:
            break
        for start, end in zip(input_poly, input_poly[1:] + input_poly[:1]):
            def side(point):
                return ((b[0] - a[0]) * (point[1] - a[1])
                        - (b[1] - a[1]) * (point[0] - a[0]))
            ds, de = side(start), side(end)
            start_in, end_in = ds >= -TOL, de >= -TOL
            if start_in:
                output.append(start)
            if start_in != end_in:
                fraction = ds / (ds - de)
                output.append((start[0] + fraction * (end[0] - start[0]),
                               start[1] + fraction * (end[1] - start[1])))
    result = []
    for point in output:
        if not result or math.dist(point, result[-1]) > TOL:
            result.append(point)
    if len(result) > 1 and math.dist(result[0], result[-1]) <= TOL:
        result.pop()
    return result


def register_mesh(cells):
    point_ids = {}
    points = []
    face_occurrences = defaultdict(list)
    materials = []
    cell_centres = []
    for cell_id, cell in enumerate(cells):
        volume = signed_volume(cell)
        if volume < -1.0e-9:
            cell = [list(reversed(face)) for face in cell]
            cells[cell_id] = cell
            volume = -volume
        if volume <= 1.0e-9:
            raise RuntimeError(f'Cell {cell_id} has non-positive volume {volume}')
        materials.append(material_id(cell))
        unique_cell_points = {tuple(round(value / TOL) for value in point): point
                              for face in cell for point in face}
        cell_centres.append(tuple(
            sum(point[i] for point in unique_cell_points.values()) / len(unique_cell_points)
            for i in range(3)))
        for face in cell:
            ids = []
            for point in face:
                key = tuple(round(v / TOL) for v in point)
                if key not in point_ids:
                    point_ids[key] = len(points)
                    points.append(point)
                ids.append(point_ids[key])
            key = tuple(sorted(ids))
            face_occurrences[key].append((cell_id, ids))
    internal = []
    unmatched = []
    for occurrences in face_occurrences.values():
        if len(occurrences) == 1:
            unmatched.append(occurrences[0])
        elif len(occurrences) == 2:
            first, second = occurrences
            internal.append((first[0], second[0], first[1]))
        else:
            raise RuntimeError(f'Non-manifold face shared by {len(occurrences)} cells')
    boundary = []
    cut_groups = defaultdict(list)
    cut_values = set(MATERIAL_PLANES)
    for occurrence in unmatched:
        cell_id, ids = occurrence
        match = None
        for axis, value in cut_values:
            if all(abs(points[point_id][axis] - value) <= TOL for point_id in ids):
                match = (axis, value)
                break
        if match is None:
            boundary.append(occurrence)
        else:
            axis, value = match
            a, b, c = (points[ids[i]] for i in range(3))
            ab = tuple(b[i] - a[i] for i in range(3))
            ac = tuple(c[i] - a[i] for i in range(3))
            normal_component = (ab[1] * ac[2] - ab[2] * ac[1],
                                ab[2] * ac[0] - ab[0] * ac[2],
                                ab[0] * ac[1] - ab[1] * ac[0])[axis]
            cut_groups[match].append((cell_id, ids, normal_component))
    for (axis, value), occurrences in cut_groups.items():
        basis = ((1, 2), (2, 0), (0, 1))[axis]
        positives = [item for item in occurrences if item[2] > 0.0]
        negatives = [item for item in occurrences if item[2] < 0.0]
        covered = defaultdict(float)
        areas = {}
        projected = {}
        for cell_id, ids, sign in occurrences:
            poly = [(points[v][basis[0]], points[v][basis[1]]) for v in ids]
            projected[(cell_id, tuple(ids))] = poly
            areas[(cell_id, tuple(ids))] = abs(polygon_area_2d(poly))
        for pos_cell, pos_ids, _ in positives:
            pos_poly = projected[(pos_cell, tuple(pos_ids))]
            for neg_cell, neg_ids, _ in negatives:
                neg_poly = projected[(neg_cell, tuple(neg_ids))]
                overlap = intersect_convex_2d(pos_poly, neg_poly)
                area = abs(polygon_area_2d(overlap)) if len(overlap) >= 3 else 0.0
                if area <= TOL * TOL:
                    continue
                face_ids = []
                for uv in overlap:
                    coordinate = [0.0, 0.0, 0.0]
                    coordinate[axis] = value
                    coordinate[basis[0]], coordinate[basis[1]] = uv
                    coordinate = tuple(coordinate)
                    key = tuple(round(v / TOL) for v in coordinate)
                    if key not in point_ids:
                        point_ids[key] = len(points)
                        points.append(coordinate)
                    face_ids.append(point_ids[key])
                if polygon_area_2d(overlap) < 0.0:
                    face_ids.reverse()
                internal.append((pos_cell, neg_cell, face_ids))
                covered[(pos_cell, tuple(pos_ids))] += area
                covered[(neg_cell, tuple(neg_ids))] += area
        for key, area in areas.items():
            if abs(covered[key] - area) > max(1.0e-5, 1.0e-5 * area):
                raise RuntimeError(f'Unmatched cut-face area on plane {axis}={value}: '
                                   f'{covered[key]} versus {area}')
    return points, internal, boundary, materials, cell_centres


HEADER = '''/* Generated from cfMesh pMesh output for Kobayashi Problem 3. */
FoamFile
{{
    version 2.0;
    format ascii;
    arch "LSB;label=32;scalar=64";
    class {class_name};
    location "constant/polyMesh";
    object {object_name};
}}
'''


def write_list(path, class_name, object_name, lines):
    with path.open('w') as stream:
        stream.write(HEADER.format(class_name=class_name, object_name=object_name))
        stream.write(f'\n{len(lines)}\n(\n')
        for line in lines:
            stream.write(f'{line}\n')
        stream.write(')\n')


def boundary_name(face, points):
    centre = tuple(sum(points[v][i] for v in face) / len(face) for i in range(3))
    distances = ((abs(centre[0]), 'xmin'), (abs(centre[0] - 60.0), 'xmax'),
                 (abs(centre[1]), 'ymin'), (abs(centre[1] - 100.0), 'ymax'),
                 (abs(centre[2]), 'zmin'), (abs(centre[2] - 60.0), 'zmax'))
    distance, name = min(distances)
    if distance > 1.0e-4:
        raise RuntimeError(f'Boundary face is not on the benchmark box: {centre}')
    return name


def write_mesh(output, points, internal, boundary, materials, cell_centres):
    output.mkdir(parents=True, exist_ok=True)
    patches = {name: [] for name in ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax')}
    bad = []
    for owner, face in boundary:
        try:
            name = boundary_name(face, points)
            axis = {'x': 0, 'y': 1, 'z': 2}[name[0]]
            expected_positive = name.endswith('max')
            a, b, c = (points[face[i]] for i in range(3))
            ab = tuple(b[i] - a[i] for i in range(3))
            ac = tuple(c[i] - a[i] for i in range(3))
            component = (ab[1] * ac[2] - ab[2] * ac[1],
                         ab[2] * ac[0] - ab[0] * ac[2],
                         ab[0] * ac[1] - ab[1] * ac[0])[axis]
            if (component > 0.0) != expected_positive:
                face = list(reversed(face))
            patches[name].append((owner, face))
        except RuntimeError:
            bad.append((owner, face))
    if bad:
        samples = [tuple(sum(points[v][i] for v in face) / len(face) for i in range(3))
                   for _, face in bad[:10]]
        raise RuntimeError(f'{len(bad)} unmatched internal faces; samples={samples}')
    oriented_internal = []
    for owner, neighbour, face in internal:
        a, b, c = (points[face[i]] for i in range(3))
        ab = tuple(b[i] - a[i] for i in range(3))
        ac = tuple(c[i] - a[i] for i in range(3))
        normal = (ab[1] * ac[2] - ab[2] * ac[1],
                  ab[2] * ac[0] - ab[0] * ac[2],
                  ab[0] * ac[1] - ab[1] * ac[0])
        outward = tuple(a[i] - cell_centres[owner][i] for i in range(3))
        if sum(normal[i] * outward[i] for i in range(3)) < 0.0:
            face = list(reversed(face))
        oriented_internal.append((owner, neighbour, face))
    internal = oriented_internal
    ordered_boundary = [item for name in patches for item in patches[name]]
    faces = [face for _, _, face in internal] + [face for _, face in ordered_boundary]
    owners = [owner for owner, _, _ in internal] + [owner for owner, _ in ordered_boundary]
    neighbours = [neighbour for _, neighbour, _ in internal]
    write_list(output / 'points', 'vectorField', 'points',
               [f'({p[0]:.16g} {p[1]:.16g} {p[2]:.16g})' for p in points])
    write_list(output / 'faces', 'faceList', 'faces',
               [f'{len(face)}(' + ' '.join(map(str, face)) + ')' for face in faces])
    write_list(output / 'owner', 'labelList', 'owner', list(map(str, owners)))
    write_list(output / 'neighbour', 'labelList', 'neighbour', list(map(str, neighbours)))
    with (output / 'boundary').open('w') as stream:
        stream.write(HEADER.format(class_name='polyBoundaryMesh', object_name='boundary'))
        stream.write('\n6\n(\n')
        start = len(internal)
        for name, values in patches.items():
            stream.write(f'{name}\n{{\n    type patch;\n    nFaces {len(values)};\n'
                         f'    startFace {start};\n}}\n')
            start += len(values)
        stream.write(')\n')
    zones = [('unused', [0])]
    zones.extend((name, [i for i, value in enumerate(materials) if value == mat_id])
                 for mat_id, name in ((1, 'source'), (2, 'duct'), (3, 'shield')))
    with (output / 'cellZones').open('w') as stream:
        stream.write(HEADER.format(class_name='regIOobject', object_name='cellZones'))
        stream.write(f'\n{len(zones)}\n(\n')
        for name, ids in zones:
            stream.write(f'{name}\n{{\n    type cellZone;\n    cellLabels List<label>\n'
                         f'    {len(ids)}\n    (\n')
            stream.writelines(f'    {cell_id}\n' for cell_id in ids)
            stream.write('    );\n}\n')
        stream.write(')\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    input_cells = load_mesh(args.input)
    cells = []
    for input_cell in input_cells:
        parts = [input_cell]
        for axis, value, transverse_ranges in CUTS:
            parts = [piece for part in parts
                     for piece in split_cell(part, axis, value, transverse_ranges)]
        cells.extend(parts)

    def is_inside(cell):
        vertices = {tuple(round(value / TOL) for value in point): point
                    for face in cell for point in face}
        centre = tuple(sum(point[i] for point in vertices.values()) / len(vertices)
                       for i in range(3))
        return (TOL < centre[0] < 60.0 - TOL
                and TOL < centre[1] < 100.0 - TOL
                and TOL < centre[2] < 60.0 - TOL)
    cells = [cell for cell in cells if is_inside(cell)]
    points, internal, boundary, materials, cell_centres = register_mesh(cells)
    write_mesh(args.output, points, internal, boundary, materials, cell_centres)
    counts = {mid: materials.count(mid) for mid in (1, 2, 3)}
    print(f'input_cells={len(input_cells)} output_cells={len(cells)} points={len(points)} '
          f'internal_faces={len(internal)} boundary_faces={len(boundary)} materials={counts}')


if __name__ == '__main__':
    main()
