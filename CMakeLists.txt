cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(FPMC)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -Wno-deprecated-declarations -pedantic -std=c++11 -g")
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/cmake) 

set(MODULES Fpmc Herwig)

#----- enable Fortran

enable_language(Fortran)

#----- external dependencies

include(UseEnvironment)
include_directories(${CLHEP_INCLUDE})

#----- build all the intermediate objects

foreach(_mod ${MODULES})
  include_directories(${PROJECT_SOURCE_DIR}/${_mod})
  add_subdirectory(${_mod})
endforeach()

#----- all useful libraries

set(LIBRARIES FPMC FPMCDeps HerwigCore ${GSL_LIB} ${GSL_CBLAS_LIB} ${LHAPDF_LIB} ${CLHEP_LIB})

#----- test executable

if(CERNLIB_LIB)
  include_directories(Examples)
  add_executable(fpmc-paw Examples/fpmc_main.f Examples/fpmc_lhe.f Examples/ntuple.f Examples/reco_ok/reco_ok.f Examples/reco_ok/pxcone_mod_new.f Examples/ffcard.f Examples/fpmc_welcome.f)
  target_link_libraries(fpmc-paw ${LIBRARIES} ${CERNLIB_LIB} ${CERNLIB_MATHLIB_LIB} ${CERNLIB_PAWLIB_LIB} ${CERNLIB_PACKLIB_LIB})
endif()

if(HEPMC_LIB)
  include_directories(${HEPMC_INCLUDE})
  message(STATUS "HepMC found at " ${HEPMC_LIB})
  include_directories(Examples)
  add_executable(fpmc-hepmc Examples/fpmc_bare.cc Examples/dummy_hwaend.f)
  target_link_libraries(fpmc-hepmc ${LIBRARIES} ${HEPMC_LIB} ${HEPMC_FIO_LIB})
endif()

#----- copies

file(COPY ${PROJECT_SOURCE_DIR}/Fpmc/External/CHIDe/Data/ DESTINATION Fpmc/External/CHIDe/Data/)
