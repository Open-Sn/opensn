sub_obj = {
  num_groups = 2,
}

#Optional parameter "limiter_type". Should create a deprecation warning
unit_testsB.TestObject.Create({
  solver_type = "B",
  coupled_field = "T",
  sub_obj1 = sub_obj,
  limiter_type = 2,
})

#Optional parameter "scheme". Should create a deprecation error.
unit_testsB.TestObject.Create({
  solver_type = "B",
  coupled_field = "T",
  sub_obj1 = sub_obj,
  scheme = "Snotty",
})
