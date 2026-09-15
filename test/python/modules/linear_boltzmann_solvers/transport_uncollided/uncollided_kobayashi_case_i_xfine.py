#!/usr/bin/env python3
"""Long Kobayashi case-i collided solve on the extra-fine (111,272-cell) tetrahedral mesh."""

from uncollided_kobayashi_utils import print_point_values, rank, run_kobayashi_case


if __name__ == "__main__":
    result = run_kobayashi_case(
        "kobayashi_dog_leg_xfine.msh",
        "i",
        "uncollided_kobayashi_case_i_xfine.h5",
        "uncollided_kobayashi_case_i_xfine.csv",
        require_existing_uncollided=True,
    )

    if rank == 0:
        print_point_values("KobayashiXfineCaseI", result["rows"])
        metrics = result["metrics"]
        print(f"KobayashiXfineExteriorMeanRatio={metrics['exterior_mean']:.8e}")
        print(f"KobayashiXfineExteriorMinRatio={metrics['exterior_min']:.8e}")
        print(f"KobayashiXfineExteriorMaxRatio={metrics['exterior_max']:.8e}")
        print(f"KobayashiXfine3CMeanRatio={metrics['mean_3c']:.8e}")

        if abs(metrics["exterior_mean"] - 1.0) > 0.05:
            raise RuntimeError(
                "Kobayashi extra-fine exterior mean reference ratio failed: "
                f"{metrics['exterior_mean']}"
            )
        if abs(metrics["mean_3c"] - 1.0) > 0.05:
            raise RuntimeError(
                f"Kobayashi extra-fine line-3C mean reference ratio failed: {metrics['mean_3c']}"
            )
        if metrics["exterior_min"] < 0.95 or metrics["exterior_max"] > 1.05:
            raise RuntimeError(
                "Kobayashi extra-fine exterior point ratios fell outside [0.95, 1.05]: "
                f"[{metrics['exterior_min']}, {metrics['exterior_max']}]"
            )
