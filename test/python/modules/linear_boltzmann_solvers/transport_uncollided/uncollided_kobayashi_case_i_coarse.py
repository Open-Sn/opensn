#!/usr/bin/env python3
"""Long Kobayashi case-i regression on the coarse tetrahedral mesh."""

from uncollided_kobayashi_utils import rank, run_kobayashi_case


if __name__ == "__main__":
    metrics = run_kobayashi_case(
        "kobayashi_dog_leg.msh",
        "i",
        "uncollided_kobayashi_coarse.h5",
        "uncollided_kobayashi_case_i_coarse.csv",
        require_existing_uncollided=True,
    )

    if rank == 0:
        print(f"KobayashiCoarseMeanRatio={metrics['mean']:.8e}")
        print(f"KobayashiCoarseMinRatio={metrics['min']:.8e}")
        print(f"KobayashiCoarseMaxRatio={metrics['max']:.8e}")
        print(f"KobayashiCoarse3CMeanRatio={metrics['mean_3c']:.8e}")
