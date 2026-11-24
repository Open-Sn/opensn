#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The console's MPIAllReduce on scalars and on lists of floats (3 ranks).

Rank r contributes the list [r, 1, -r, 0.5 r]. The sums, maxima and minima over n ranks follow in
closed form, so each result is checked exactly. The list form needs pybind11's STL conversions,
without which every list call raises TypeError. Under the Python module, mpi4py stands in for the
console function, so the test checks only the closed forms there.
"""

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        if isinstance(value, list):
            return [MPI.COMM_WORLD.allreduce(v, op=ops[op]) for v in value]
        return MPI.COMM_WORLD.allreduce(value, op=ops[op])

if __name__ == "__main__":
    n = size
    local = [float(rank), 1.0, -float(rank), 0.5 * rank]
    rank_sum = n * (n - 1) / 2
    expected = {
        "sum": [rank_sum, float(n), -rank_sum, 0.5 * rank_sum],
        "max": [n - 1.0, 1.0, 0.0, 0.5 * (n - 1)],
        "min": [0.0, 1.0, -(n - 1.0), 0.0],
    }
    passed = 0
    for op, values in expected.items():
        result = MPIAllReduce(local, op)
        passed += int(isinstance(result, list) and list(result) == values)
    passed += int(MPIAllReduce(rank, "sum") == int(rank_sum))
    passed += int(MPIAllReduce(0.5 * rank, "max") == 0.5 * (n - 1))
    passed += int(MPIAllReduce([], "sum") == [])
    if rank == 0:
        print(f"CONSOLE_ALLREDUCE_PASSED= {passed}")
