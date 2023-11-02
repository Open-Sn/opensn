sub_obj = {
  num_groups = 2,
}

chi_unit_testsB.ChildTestObject.Create(
  {
    solver_type = "C",
    coupled_field = "T",
    sub_obj1 = sub_obj,
    num_sub_groups = 3
  })

--Required parameter "format". Should create a deprecation error.
chi_unit_testsB.TestObject.Create(
  {
    solver_type = "B",
    coupled_field = "T",
    sub_obj1 = sub_obj,
    format = true
  })