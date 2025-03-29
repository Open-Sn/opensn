if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.xs import MultiGroupXS

# Create cross sections
xss = []
for m in range(7):
    xss.append( MultiGroupXS() )

xss[0].LoadFromOpenSn("c5g7/materials/XS_water.xs")
xss[1].LoadFromOpenSn("c5g7/materials/XS_UO2.xs")
xss[2].LoadFromOpenSn("c5g7/materials/XS_7pMOX.xs")
xss[3].LoadFromOpenSn("c5g7/materials/XS_guide_tube.xs")
xss[4].LoadFromOpenSn("c5g7/materials/XS_4_3pMOX.xs")
xss[5].LoadFromOpenSn("c5g7/materials/XS_8_7pMOX.xs")
xss[6].LoadFromOpenSn("c5g7/materials/XS_fission_chamber.xs")

num_groups = xss[0].num_groups
print("Num groups: ", num_groups)

# Create materials
xs_map = []
for m in range(0, len(xss)):
    xs_map.append( { "block_ids": [ m ], "xs": xss[m] } )
