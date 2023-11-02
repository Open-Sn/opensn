include_directories("LinearBoltzmannSolvers")

add_subdirectory("DiffusionSolver")
add_subdirectory("LinearBoltzmannSolvers")
add_subdirectory("CFEMDiffusion")
add_subdirectory("DFEMDiffusion")
add_subdirectory("MGDiffusion")
add_subdirectory("FVDiffusion")
add_subdirectory("PointReactorKinetics")

set(MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/lua/chi_modules_lua.cc")

set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)
