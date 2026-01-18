#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# --- Lebedev Quadrature Test 1 ---
# Create a Lebedev quadrature of order 7
print("\n--- Testing Lebedev Quadrature, Order 7 ---")
quad1_sum = sum(LebedevQuadrature3DXYZ(quadrature_order=7, scattering_order=0).weights)
print(f"Weight-Sum-1={quad1_sum:.3e}\n\n")

# --- Lebedev Quadrature Test 2 ---
# Create a higher order Lebedev quadrature
print("\n--- Testing Lebedev Quadrature, Order 15 ---")
quad2_sum = sum(LebedevQuadrature3DXYZ(quadrature_order=15, scattering_order=0).weights)
print(f"Weight-Sum-2={quad2_sum:.3e}\n\n")
