"""Shared data for OECD/NEA Kobayashi benchmark Problem 3."""


MESH_FILENAME = "kobayashi_dog_leg.msh"
UNCOLLIDED_FILENAME = "kobayashi_dog_leg_uncollided.h5"
RESULTS_FILENAME = "kobayashi_dog_leg_results.csv"

MATERIAL_IDS = {
    "source": 1,
    "duct": 2,
    "shield": 3,
}

SOURCE_STRENGTH = 1.0
SOURCE_VOLUME = 10.0**3

# Problem 3 reference scalar fluxes from OECD/NEA NSC/DOC(2000)4.
REFERENCE_POINTS = {
    "3A": [
        ((5.0, 5.0, 5.0), 5.95659e0, 8.61578e0),
        ((5.0, 15.0, 5.0), 1.37185e0, 2.16130e0),
        ((5.0, 25.0, 5.0), 5.00871e-1, 8.93784e-1),
        ((5.0, 35.0, 5.0), 2.52429e-1, 4.78052e-1),
        ((5.0, 45.0, 5.0), 1.50260e-1, 2.89424e-1),
        ((5.0, 55.0, 5.0), 9.91726e-2, 1.92698e-1),
        ((5.0, 65.0, 5.0), 4.22623e-2, 1.04982e-1),
        ((5.0, 75.0, 5.0), 1.14703e-2, 3.37544e-2),
        ((5.0, 85.0, 5.0), 3.24662e-3, 1.08158e-2),
        ((5.0, 95.0, 5.0), 9.48324e-4, 3.39632e-3),
    ],
    "3B": [
        ((5.0, 55.0, 5.0), 9.91726e-2, 1.92698e-1),
        ((15.0, 55.0, 5.0), 2.45041e-2, 6.72147e-2),
        ((25.0, 55.0, 5.0), 4.54477e-3, 2.21799e-2),
        ((35.0, 55.0, 5.0), 1.42960e-3, 9.90646e-3),
        ((45.0, 55.0, 5.0), 2.64846e-4, 3.39066e-3),
        ((55.0, 55.0, 5.0), 9.14210e-5, 1.05629e-3),
    ],
    "3C": [
        ((5.0, 95.0, 35.0), 3.27058e-5, 3.44804e-4),
        ((15.0, 95.0, 35.0), 2.68415e-5, 2.91825e-4),
        ((25.0, 95.0, 35.0), 1.70019e-5, 2.05793e-4),
        ((35.0, 95.0, 35.0), 3.37981e-5, 2.62086e-4),
        ((45.0, 95.0, 35.0), 6.04893e-6, 1.05367e-4),
        ((55.0, 95.0, 35.0), 3.36460e-6, 4.44962e-5),
    ],
}


def cross_sections(case):
    """Return ``(sigma_t, scattering_ratio)`` by material for case i or ii."""
    if case == "i":
        return {
            "source": (0.1, 0.0),
            "duct": (1.0e-4, 0.0),
            "shield": (0.1, 0.0),
        }
    if case == "ii":
        return {
            "source": (0.1, 0.5),
            "duct": (1.0e-4, 0.5),
            "shield": (0.1, 0.5),
        }
    raise ValueError(f"KOBAYASHI_CASE must be 'i' or 'ii', got {case!r}")
