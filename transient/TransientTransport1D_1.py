#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D transient k-eigenvalue problem
"""

import os
import sys
import math

if __name__ == "__main__":

    nodes = []
    N = 100
    L = 30.0
    dx = L / N
    nodes = [i * dx for i in range(N + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    num_groups = 1
    test_material = MultiGroupXS();
    test_material.LoadFromOpenSn("xs_inf_critical_1g.cxs")
    #test_material2 = MultiGroupXS();
    #test_material.LoadFromOpenSn("subcritical_1g.cxs")
    #materials = [ ( test_material1, 0.002644086 ), ( test_material2, 0.0424459 ) ]
    #test_material = MultiGroupXS()
    #test_material.Combine(materials)
    #xs = chiPhysicsTransportXSMakeCombined({{xs, 0.00264086}}) -- just sub-critical
    #xs = chiPhysicsTransportXSMakeCombined({{xs, 0.0424459}}) -- just sub-critical

    pquad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=0)

    # Create solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups-1),
                "angular_quadrature": pquad,
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {
                "block_ids": [0],
                "xs": test_material
            }
        ],
        options={
            "use_precursors": True
        },
        time_dependent=True
    )
    solver = TransientKEigenSolver(problem=phys)
    solver.Initialize()
    solver.Execute()

#chiLBTSSetProperty(phys1, "TIMESTEP", 1e-1)
#chiLBTSSetProperty(phys1, "VERBOSITY_LEVEL", 0)
#chiLBTSSetProperty(phys1, "TIMESTEP_METHOD", "CRANK_NICHOLSON")

#time = 0.0
#time_stop = 20.0
#k=0
#while (time < time_stop) do
#    k = k + 1
#    solver.Step(phys1)
#    FRf = ComputeFissionRate(phys1,"NEW")
#    FRi = ComputeFissionRate(phys1,"OLD")
#    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
#    time = chiLBTSGetProperty(phys1, "TIME")
#    period = dt/math.log(FRf/FRi)
#end
