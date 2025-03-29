if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.xs import MultiGroupXS

xss = []
xss.append( MultiGroupXS() )
xss[0].LoadFromOpenSn("../transport_keigen/xs_water_g2.xs")
xss.append( MultiGroupXS() )
xss[1].LoadFromOpenSn("../transport_keigen/xs_fuel_g2.xs")
num_groups = xss[0].num_groups