#!/usr/bin/env python3
"""Generate the extra-fine Kobayashi Problem 3 uncollided file for long regressions."""

from uncollided_kobayashi_utils import generate_kobayashi_uncollided, rank


if __name__ == "__main__":
    generate_kobayashi_uncollided(
        "kobayashi_dog_leg_xfine.msh",
        "uncollided_kobayashi_case_i_xfine.h5",
        case="i",
    )

    if rank == 0:
        print("KobayashiXfineUncollidedGenerated=1")
