#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script reads quadrature point files, sums a specific column of weights,
and prints normalized weight sums. It demonstrates creating a quadrature,
printing its points to a file, refining the quadrature locally, and cleaning up.
"""

import os
import math


def split_string(input_str, separator):
    """
    Splits the input_str using the specified separator.
    """
    # Using str.split() directly is sufficient for simple splits.
    return input_str.split(separator)


def get_weight_sum(quad_points_file):
    """
    Opens a quadrature points file, reads each line,
    splits the line to extract the 4th value (index 3) as a float,
    and sums these values.
    """
    print(f"\nReading: {quad_points_file}")
    weight_sum = 0.0
    try:
        with open(quad_points_file, "r") as file:
            file.readline()
            for line in file:
                # Remove any extra whitespace and split by space.
                # (Alternatively, you could simply use line.split() to split on any whitespace.)
                values = split_string(line.strip(), ",")
                weight_sum += float(values[3])
    except Exception:
        print(f"Error: Could not open file {quad_points_file}")
    print(weight_sum)
    return weight_sum


# --- Quadrature-1 ---
# Create a quadrature with initial refinement level of 0.
cquad1 = SLDFESQuadrature(level=0)
cquad1.PrintQuadratureToFile("TestQuad1")
quad1_sum = get_weight_sum("TestQuad1_points.csv")
print(f"Weight-Sum-1={(quad1_sum / (4 * math.pi)):.3e}\n\n")

# --- Quadrature-2 ---
cquad2 = SLDFESQuadrature(level=2)
# Locally refine the quadrature near a given point with a given angular spread.
cquad2.LocallyRefine(Vector3(0.25, -0.85, 1.0), 30.0 * math.pi / 180, False)
cquad2.PrintQuadratureToFile("TestQuad2")
quad2_sum = get_weight_sum("TestQuad2_points.csv")
print(f"Weight-Sum-2={(quad2_sum / (4 * math.pi)):.3e}\n\n")

# Clean up: remove generated files.
os.system("rm TestQuad1* TestQuad2*")
