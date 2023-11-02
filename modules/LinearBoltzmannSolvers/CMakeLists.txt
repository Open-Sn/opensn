add_subdirectory("A_LBSSolver")
add_subdirectory("B_DiffusionDFEMSolver")
add_subdirectory("B_DiscreteOrdinatesSolver")
add_subdirectory("C_DiscreteOrdinatesAdjointSolver")
add_subdirectory("C_DiscreteOrdinatesCurvilinearSolver")
add_subdirectory("Executors")

#add_subdirectory("Ca_DO_SteadyState")
#add_subdirectory("D_DO_Transient")



set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)