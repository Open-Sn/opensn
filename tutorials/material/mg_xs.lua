--[[ @doc
# Multigroup Cross Sections

## Background

OpenSn is not provided with cross-section libraries. Users are expected to supply their own multigroup cross-section data.
One may use other open-source software to generate this data (e.g., NJOY, Dragon, OpenMC).

The cross sections are read from a plain text file. The ```OPENSN_XSFILE``` format  of that file is as follows:
```
# Add descriptive comments as appropriate
NUM_GROUPS ng
NUM_MOMENTS nmom

SIGMA_T_BEGIN
0 value
.
.
ng-1 value
SIGMA_T_END

SIGMA_A_BEGIN
0 value
.
.
ng-1 value
SIGMA_A_END

TRANSFER_MOMENTS_BEGIN
# Add descriptive comments as appropriate
M_GPRIME_G_VAL 0 0 0 value
.
M_GPRIME_G_VAL moment gprime g value
.
M_GPRIME_G_VAL nmom-1 ng-1 ng-1 value
TRANSFER_MOMENTS_END

```

## Cross Sections

We load the cross sections from OpenSn file format.
--]]
xs_matA = xs.LoadFromOpenSn("xs_1g_MatA.xs")
