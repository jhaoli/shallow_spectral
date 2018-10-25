cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(stswm LANGUAGES Fortran)

find_package(MPI)
if (NOT MPI_Fortran_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DNO_MPI")
else ()
  include_directories(${MPI_Fortran_INCLUDE_PATH})
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -ffpe-summary=none -ffree-line-length-none -Ofast -fopenmp -fno-align-commons")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -fast -fPIC")
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


add_subdirectory(lib)
#include_directories(${CMAKE_BINARY_DIR}/fortran_container)


set(source_files
  src/advect.f
  src/agpwrt.f
  src/anlytc.f
  src/bf2.f
  src/bubfnc.f
  src/calp.f
  src/ceps.f
  src/comp1.f
  src/cpcnrc.f
  src/cpmpxy.f
  src/cpmvxy.f
  src/d2bubf.f
  src/dbubf.f
  src/dcalp.f
  src/dceps.f
  src/dft991.f
  src/dftfax.f
  src/drawcl.f
  src/dtrns.f
  src/duv.f
  src/epslon.f
  src/erranl.f
  src/eval.f
  src/force.f
  src/ftrnde.f
  src/ftrndi.f
  src/ftrnex.f
  src/ftrnim.f
  src/ftrnpe.f
  src/ftrnpi.f
  src/ftrnve.f
  src/fu.f
  src/func2.f
  src/glat.f
  src/glats.f
  src/glon.f
  src/infld.f
  src/init.f
  src/inptp.f
  src/input.f
  src/labtop.f
  src/nrgtcs.f
  src/ordleg.f
  src/ortho.f
  src/outptp.f
  src/plots.f
  src/prnt.f
  src/rotate.f
  src/shader.f
  src/shtrns.f
  src/spcanl.f
  src/step.f
  src/stswm.f
  src/us.f
  src/uv.f
  src/velvct.f
  src/weight.f
  src/zd.f
  lib/adquad.f
  lib/fft99f.f
  lib/ncarg.f
  lib/netcdf.f
)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/.mods")

add_executable(stswm_test.exe ${source_files})
target_link_libraries(stswm_test.exe
  ${MPI_Fortran_LIBRARIES} fftpack netcdf netcdff)
if (MPI_Fortran_FOUND)
  set_target_properties(stswm_test.exe PROPERTIES
    LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
  )
endif ()
