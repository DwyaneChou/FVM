cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(fvm LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  #set(CPU_ARCH HASWELL)
  set(CPU_ARCH BROADWELL)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -r8 -qopenmp -mkl -x${CPU_ARCH} -ax${CPU_ARCH} -arch ${CPU_ARCH} -mcpu=${CPU_ARCH} -m64")
  if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -check all -check noarg_temp_created")
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  if (APPLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -floop-parallelize-all -fdefault-real-8 -ffree-line-length-none")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -floop-parallelize-all -fdefault-real-8 -ftree-parallelize-loops=8 -ffree-line-length-none -lgomp -fopenmp")
  endif ()
  if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -gdwarf -Og -fbacktrace -fcheck=all")
  endif ()
endif ()

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

set(source_files
  src/fvm.F90
  src/constants_mod.F90
  src/parameters_mod.F90
  src/projection_mod.F90
  src/mesh_mod.F90
  src/stat_mod.F90
  src/tend_mod.F90
  src/math_mod.F90
  src/reconstruction_mod.F90
  src/qr_solver_mod.F90
  src/test_case_mod.F90
  src/io_mod.F90
  src/spatial_operators_mod.F90
  src/temporal_mod.F90
  src/diag_mod.F90
  src/quadpack.F90
  src/weno_coef_mod.F90
)

add_library(fvm ${source_files})
target_link_libraries(fvm netcdf netcdff)

add_executable(fvm.exe src/fvm.F90)
target_link_libraries(fvm.exe fvm)

