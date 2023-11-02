sub_obj = {
  num_groups = 2,
}

--Required parameter "use_ragusas_stuff". Should create an error with
--                                        renaming message
chi_unit_testsB.TestObject.Create(
  {
    solver_type = "B",
    coupled_field = "T",
    sub_obj1 = sub_obj,
    use_ragusas_stuff = true
  })