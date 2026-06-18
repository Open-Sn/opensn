#!/usr/bin/env python3
"""Long Kobayashi case-i regression on the fine tetrahedral mesh."""

from uncollided_kobayashi_utils import print_point_values, rank, run_kobayashi_case


if __name__ == "__main__":
    result = run_kobayashi_case(
        "kobayashi_dog_leg_fine.msh",
        "i",
        "uncollided_kobayashi_case_i_fine.h5",
        "uncollided_kobayashi_case_i_fine.csv",
        require_existing_uncollided=True,
    )

    if rank == 0:
        print_point_values("KobayashiFineCaseI", result["rows"])
        metrics = result["metrics"]
        print(f"KobayashiFineMeanRatio={metrics['mean']:.8e}")
        print(f"KobayashiFineMinRatio={metrics['min']:.8e}")
        print(f"KobayashiFineMaxRatio={metrics['max']:.8e}")
        print(f"KobayashiFine3CMeanRatio={metrics['mean_3c']:.8e}")
