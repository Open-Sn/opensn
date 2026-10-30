# One group, sigma_t = 200 /cm, scattering ratio c = 0.99, Henyey-Greenstein
# scattering with mean cosine g = 0.9 truncated at P3: sigma_s,l = c sigma_t g^l.
NUM_GROUPS 1
NUM_MOMENTS 4

SIGMA_T_BEGIN
0 200.0
SIGMA_T_END

SIGMA_A_BEGIN
0 2.0
SIGMA_A_END

TRANSFER_MOMENTS_BEGIN
M_GFROM_GTO_VAL 0 0 0 198.0
M_GFROM_GTO_VAL 1 0 0 178.2
M_GFROM_GTO_VAL 2 0 0 160.38
M_GFROM_GTO_VAL 3 0 0 144.342
TRANSFER_MOMENTS_END
