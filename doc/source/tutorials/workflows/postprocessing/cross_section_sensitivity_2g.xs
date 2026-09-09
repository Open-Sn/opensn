# Two-group downscattering cross sections for the sensitivity tutorial.
NUM_GROUPS 2
NUM_MOMENTS 1

SIGMA_T_BEGIN
0 1.0
1 0.8
SIGMA_T_END

TRANSFER_MOMENTS_BEGIN
# M_GFROM_GTO_VAL moment from_group to_group value
M_GFROM_GTO_VAL 0 0 0 0.2
M_GFROM_GTO_VAL 0 0 1 0.3
M_GFROM_GTO_VAL 0 1 0 0.0
M_GFROM_GTO_VAL 0 1 1 0.1
TRANSFER_MOMENTS_END
