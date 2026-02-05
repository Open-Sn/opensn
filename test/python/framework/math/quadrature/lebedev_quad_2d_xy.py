#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test 2D Lebedev quadrature (upper hemisphere only)

# --- Lebedev Quadrature 2D Test 1 ---
# Create a 2D Lebedev quadrature of order 7
print("\n--- Testing Lebedev Quadrature 2D, Order 7 ---")
quad1_sum = sum(LebedevQuadrature2DXY(quadrature_order=7, scattering_order=0).weights)
print(f"Weight-Sum-1={quad1_sum:.3e}\n\n")

# --- Lebedev Quadrature 2D Test 2 ---
# Create a higher order 2D Lebedev quadrature
print("\n--- Testing Lebedev Quadrature 2D, Order 15 ---")
quad2_sum = sum(LebedevQuadrature2DXY(quadrature_order=15, scattering_order=0).weights)
print(f"Weight-Sum-2={quad2_sum:.3e}\n\n")
