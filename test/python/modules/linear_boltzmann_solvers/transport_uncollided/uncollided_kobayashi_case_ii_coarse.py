#!/usr/bin/env python3
"""Long Kobayashi case-ii regression on the coarse tetrahedral mesh."""

from uncollided_kobayashi_utils import print_point_values, rank, run_kobayashi_case


if __name__ == "__main__":
    result = run_kobayashi_case(
        "kobayashi_dog_leg.msh",
        "ii",
        "uncollided_kobayashi_coarse.h5",
        "uncollided_kobayashi_case_ii_coarse.csv",
        require_existing_uncollided=True,
    )

    if rank == 0:
        print_point_values("KobayashiCoarseCaseII", result["rows"])
        metrics = result["metrics"]
        print(f"KobayashiCoarseCaseIIMeanRatio={metrics['mean']:.8e}")
        print(f"KobayashiCoarseCaseIIMinRatio={metrics['min']:.8e}")
        print(f"KobayashiCoarseCaseIIMaxRatio={metrics['max']:.8e}")
        print(f"KobayashiCoarseCaseII3CMeanRatio={metrics['mean_3c']:.8e}")

        if abs(metrics["exterior_mean"] - 1.0) > 0.08:
            raise RuntimeError(
                "Kobayashi case-II exterior mean reference ratio failed: "
                f"{metrics['exterior_mean']}"
            )
        if abs(metrics["mean_3c"] - 1.0) > 0.08:
            raise RuntimeError(
                f"Kobayashi case-II line-3C mean reference ratio failed: {metrics['mean_3c']}"
            )
