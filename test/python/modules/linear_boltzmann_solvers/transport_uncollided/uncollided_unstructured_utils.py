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


def uncollided_2d(strength, sigma_t, source, point):
    radius = math.dist(source[:2], point[:2])
    return strength * math.exp(-sigma_t * radius) / (2.0 * math.pi * radius)


def uncollided_3d(strength, sigma_t, source, point):
    radius = math.dist(source, point)
    return strength * math.exp(-sigma_t * radius) / (4.0 * math.pi * radius * radius)


def uncollided_escape_2d_rectangle(sigma_t, source, bounds, num_angles=100000):
    xmin, xmax, ymin, ymax = bounds
    escape = 0.0
    for index in range(num_angles):
        angle = 2.0 * math.pi * (index + 0.5) / num_angles
        omega_x = math.cos(angle)
        omega_y = math.sin(angle)
        distances = []
        if omega_x > 0.0:
            distances.append((xmax - source[0]) / omega_x)
        elif omega_x < 0.0:
            distances.append((xmin - source[0]) / omega_x)
        if omega_y > 0.0:
            distances.append((ymax - source[1]) / omega_y)
        elif omega_y < 0.0:
            distances.append((ymin - source[1]) / omega_y)
        distance = min(value for value in distances if value > 0.0)
        escape += math.exp(-sigma_t * distance)
    return escape / num_angles


def relative_error(value, reference):
    return abs(value - reference) / abs(reference)


def remove_file(file_name):
    try:
        os.remove(file_name)
    except FileNotFoundError:
        pass
