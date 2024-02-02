

block =
{
  enabled = true,
  it_method = "gmres",
  nl_abs_tol = 1.0e-12,
  nl_max_its = 33,
  sub1 =
  {
    ax_method = 2,
    l_abs_tol = 1.0e-2
  },
  sub2 =
  {
    ax_method = 3,
    l_abs_tol = 1.0e-3,
    blocks = {99, 98, 97},
    cblocks = {{1,2,3},{4,5,6},{7,8,9}}
  }
}

unit_tests.ParameterBlock_Test00(--[[verbose=]]true, block)
