#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D CMFD correction limiter stress test with intentionally over-relaxed corrections.

This regression reuses the standard 1D 1-group CMFD k-eigenvalue problem, but
sets an intentionally large CMFD relaxation factor and disables damping retries
by requiring a minimum damping of 1.0. The first CMFD correction should be
rejected by the correction limiter, the unaccelerated transport update should be
kept, and the verbose CMFD diagnostics should report a skipped correction
without nonfinite scalar fluxes or an invalid k-eigenvalue.
"""

import os
import runpy

cmfd_verbose = True
cmfd_limiter_stress = True
cmfd_relaxation = 100.0
cmfd_correction_max_attempts = 2
cmfd_correction_min_damping = 1.0
pi_solver_max_iters = 2
cmfd_solver_max_iters = 3

results = runpy.run_path(
    os.path.join(os.path.dirname(__file__), "keigenvalue_transport_1d_1g_cmfd.py"),
    init_globals=globals(),
    run_name="__main__",
)

assert results["pi_iterations"] == pi_solver_max_iters
assert results["cmfd_iterations"] == cmfd_solver_max_iters
