# source list
if(OPENSN_WITH_CUDA OR OPENSN_WITH_HIP)
  add_library(libcaribou INTERFACE)
  target_include_directories(
    libcaribou
    INTERFACE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:include>
  )
elseif(OPENSN_WITH_SYCL)
  file(GLOB LIBCARIBOU_SRCS CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/caribou/sycl/*.cpp
  )
  add_library(libcaribou SHARED ${LIBCARIBOU_SRCS})
  set_target_properties(libcaribou PROPERTIES
    OUTPUT_NAME crb
  )
  target_include_directories(
    libcaribou
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
  )
  target_compile_features(libcaribou PUBLIC cxx_std_20)
  target_compile_options(libcaribou PRIVATE -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend "-device pvc")
  target_link_options(libcaribou PRIVATE -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend "-device pvc")
endif()
