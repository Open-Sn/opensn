#!/usr/bin/env python3

import math
import os


def mesh_path(file_name):
    return os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../../../assets/mesh", file_name)
    )


def point_value(field_function, point, interpolation_type, vector_type):
    interpolation = interpolation_type()
    interpolation.SetPointOfInterest(vector_type(*point))
    interpolation.AddFieldFunction(field_function)
    interpolation.Execute()
    return interpolation.GetPointValue()


def volume_integral(field_function, logical_volume, interpolation_type):
    interpolation = interpolation_type()
    interpolation.SetOperationType("sum")
    interpolation.SetLogicalVolume(logical_volume)
    interpolation.AddFieldFunction(field_function)
    interpolation.Execute()
    return interpolation.GetValue()


def volume_minimum(field_function, logical_volume, interpolation_type):
    interpolation = interpolation_type()
    interpolation.SetOperationType("min")
    interpolation.SetLogicalVolume(logical_volume)
    interpolation.AddFieldFunction(field_function)
    interpolation.Execute()
    return interpolation.GetValue()


def uncollided_3d(strength, sigma_t, source, point):
    radius = math.dist(source, point)
    return strength * math.exp(-sigma_t * radius) / (4.0 * math.pi * radius * radius)


def extruded_mesh(extruder_type, file_mesh_type, file_name, height=2.0, num_layers=8):
    """Prism mesh from extruding a 2D mesh file over z in [0, height]."""
    return extruder_type(
        inputs=[file_mesh_type(filename=mesh_path(file_name))],
        layers=[{"z": height, "n": num_layers}],
    ).Execute()


def gauss_legendre(order):
    """Gauss-Legendre points and weights on [-1, 1] by Newton iteration."""
    points, weights = [], []
    for i in range(order):
        x = math.cos(math.pi * (i + 0.75) / (order + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, order + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = order * (x * p1 - p0) / (x * x - 1.0)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1.0e-15:
                break
        points.append(x)
        weights.append(2.0 / ((1.0 - x * x) * dp * dp))
    return points, weights


GL_POINTS, GL_WEIGHTS = gauss_legendre(16)


def composite_rule(a, b, num_panels):
    width = (b - a) / num_panels
    rule = []
    for p in range(num_panels):
        left = a + p * width
        for x, w in zip(GL_POINTS, GL_WEIGHTS):
            rule.append((left + 0.5 * width * (x + 1.0), 0.5 * width * w))
    return rule


def box_angular_integral(integrand, source, bounds, num_panels=20):
    """Integral over all directions of integrand(R) for a source in an axis-aligned box.

    R is the distance from the source to the box boundary along the direction. Each face is
    integrated in Cartesian coordinates with dOmega = h / r^3 dA, where h is the distance from
    the source to the face plane. Faces whose plane contains the source subtend no solid angle.
    """
    total = 0.0
    for axis in range(3):
        for side in range(2):
            h = abs(bounds[axis][side] - source[axis])
            if h == 0.0:
                continue
            a, b = [t for t in range(3) if t != axis]
            rule_a = composite_rule(bounds[a][0] - source[a], bounds[a][1] - source[a], num_panels)
            rule_b = composite_rule(bounds[b][0] - source[b], bounds[b][1] - source[b], num_panels)
            for x, wx in rule_a:
                for y, wy in rule_b:
                    r = math.sqrt(h * h + x * x + y * y)
                    total += wx * wy * h / r**3 * integrand(r)
    return total


def uncollided_escape_3d_box(sigma_t, source, bounds):
    """Fraction of a point source's emission that leaves a homogeneous box uncollided."""
    return box_angular_integral(lambda r: math.exp(-sigma_t * r), source, bounds) / (4.0 * math.pi)


def single_scatter_point_response(sigma_s, sigma_t_from, sigma_t_to, distance):
    """Infinite-medium scalar flux at B in group h from a unit isotropic point source at A.

    The source emits in group g, a particle reaches group h only by one isotropic g -> h
    scattering with cross section sigma_s, and neither group scatters otherwise. With
    r1 = |r - A|, r2 = |r - B|, and d = |A - B|, the flux is the double integral

        sigma_s / (8 pi d) int int exp(-sigma_t_from r1 - sigma_t_to r2) / (r1 r2) dr1 dr2

    over |r1 - r2| <= d <= r1 + r2 (prolate spheroidal coordinates, integrated over the
    azimuth). It is symmetric in the two total cross sections, as reciprocity requires. With
    u = r1 + r2 and v = r1 - r2, 1 / (u^2 - v^2) splits into 1 / (u - v) + 1 / (u + v), and
    the substitutions log(u -+ v) and log(u - d) remove the endpoint singularities.
    """
    a = 0.5 * (sigma_t_from + sigma_t_to)
    b = 0.5 * (sigma_t_from - sigma_t_to)
    d = distance

    def inner(u):
        total = 0.0
        rule = composite_rule(math.log(u - d), math.log(u + d), 4)
        for sign in (1.0, -1.0):
            for t, w in rule:
                total += w * math.exp(-b * sign * (u - math.exp(t)))
        return math.exp(-a * u) * total / u

    total = 0.0
    for tau, w in composite_rule(-30.0, math.log(60.0 / a), 40):
        xi = math.exp(tau)
        total += w * xi * inner(d + xi)
    return sigma_s * total / (8.0 * math.pi * d)


def relative_error(value, reference):
    return abs(value - reference) / abs(reference)


def remove_file(file_name):
    try:
        os.remove(file_name)
    except FileNotFoundError:
        pass
