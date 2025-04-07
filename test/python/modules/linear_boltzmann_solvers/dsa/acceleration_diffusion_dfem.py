#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import math

# Setup mesh
nodes = []
N = 10
L = 2.0
xmin = -L / 2
dx = L / N
for i in range(N + 1):
    nodes.append(xmin + i * dx)
meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
grid = meshgen.Execute()
grid.SetUniformBlockID(0)
grid.SetOrthogonalBoundaries()


def MMS_phi(pt):
    return math.cos(math.pi * pt.x) + math.cos(math.pi * pt.y)


def MMS_q(pt):
    return math.pi**2 * (math.cos(math.pi * pt.x) + math.cos(math.pi * pt.y))


mms_phi_fn = ScalarSpatialFunction(MMS_phi)
mms_q_fn = ScalarSpatialFunction(MMS_q)
acceleration_Diffusion_DFEM(grid)
MPIBarrier()
if rank == 0:
    os.system("rm SimTest_92*")
