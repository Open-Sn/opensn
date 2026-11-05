add_library(opensn_precice SHARED
  cell_power.cc
  precice_adapter.cc
)

target_include_directories(opensn_precice
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/opensn>
)

target_link_libraries(opensn_precice
  PUBLIC
    libopensn
  PRIVATE
    precice::precice
)

target_compile_options(opensn_precice PRIVATE ${OPENSN_CXX_FLAGS})
set_target_properties(opensn_precice PROPERTIES
  EXPORT_NAME precice
  OUTPUT_NAME opensnprecice
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

install(
  TARGETS opensn_precice
  EXPORT openSnTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
