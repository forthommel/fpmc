#----- external dependencies

if($ENV{HOSTNAME} MATCHES "^lxplus[0-9]+.cern.ch")
  set(BASE_DIR "/cvmfs/sft.cern.ch/lcg/external")
  set(CMAKE_CXX_COMPILER "${BASE_DIR}/gcc/6.1.0/x86_64-slc6/bin/g++")
  set(CMAKE_Fortran_COMPILER "${BASE_DIR}/gcc/6.1.0/x86_64-slc6/bin/gfortran")

  set(CERNLIB_DIR "${BASE_DIR}/cernlib/2006a/x86_64-slc6-gcc47-opt")
  set(GSL_DIR "${BASE_DIR}/GSL/1.14/x86_64-slc5-gcc44-opt")
  set(CLHEP_DIR "${BASE_DIR}/clhep/2.2.0.4/x86_64-slc6-gcc48-opt")
  set(HEPMC_DIR "${BASE_DIR}/HepMC/2.06.08/x86_64-slc6-gcc48-opt")
  set(LHAPDF_DIR "${BASE_DIR}/MCGenerators/lhapdf/5.8.9/x86_64-slc6-gcc46-opt")
  set(PYTHIA8_DIR "${BASE_DIR}/MCGenerators_lcgcmt67c/pythia8/230/x86_64-slc6-gcc48-opt")

  message(STATUS "Compiling on LXPLUS. Do not forget to source the environment variables!")
  message(STATUS "    e.g. `source /cvmfs/sft.cern.ch/lcg/external/gcc/6.1.0/x86_64-slc6/setup.sh`")

  #find_library(CERNLIB_LIB NAMES cernlib HINTS "${CERNLIB_DIR}/lib")
  find_library(CERNLIB_MATHLIB_LIB NAMES mathlib HINTS "${CERNLIB_DIR}/lib")
  find_library(CERNLIB_PACKLIB_LIB NAMES packlib HINTS "${CERNLIB_DIR}/lib")
  find_library(CERNLIB_PAWLIB_LIB NAMES pawlib HINTS "${CERNLIB_DIR}/lib")
  if(CERNLIB_LIB)
    message(STATUS "CERNLIB found in ${CERNLIB_LIB}")
  endif()
  if(CERNLIB_PAWLIB_LIB)
    message(STATUS "PAW found in ${CERNLIB_PAWLIB_LIB}")
  endif()

  find_library(GSL_LIB gsl HINTS "${GSL_DIR}/lib")
  find_library(GSL_CBLAS_LIB gslcblas HINTS "${GSL_DIR}/lib")
  #--- searching for CLHEP
  find_library(CLHEP_LIB CLHEP HINTS "${CLHEP_DIR}/lib")
  find_path(CLHEP_INCLUDE CLHEP HINTS "${CLHEP_DIR}/include")
  #--- searching for LHAPDF
  find_library(LHAPDF_LIB LHAPDF HINTS "${LHAPDF_DIR}/lib")
  #--- searching for HepMC
  find_library(HEPMC_LIB HepMC HINTS "${HEPMC_DIR}/lib")
  find_library(HEPMC_FIO_LIB HepMCfio HINTS "${HEPMC_DIR}/lib")
  find_path(HEPMC_INCLUDE HepMC HINTS "${HEPMC_DIR}/include")
else()
  find_library(GSL_LIB gsl)
  find_library(GSL_CBLAS_LIB gslcblas)
  find_library(CLHEP_LIB CLHEP)
  find_path(CLHEP_INCLUDE CLHEP)
  find_library(LHAPDF_LIB LHAPDF)
  find_library(HEPMC_LIB HepMC)
  find_library(HEPMC_FIO_LIB HepMCfio)
  find_path(HEPMC_INCLUDE HepMC)
endif()
#--- searching for Pythia 8
set(PYTHIA8_DIRS $ENV{PYTHIA8_DIR} ${PYTHIA8_DIR} /usr /usr/local /opt/pythia8)
find_library(PYTHIA8 pythia8 HINTS ${PYTHIA8_DIRS} PATH_SUFFIXES lib)
find_path(PYTHIA8_INCLUDE Pythia8 HINTS ${PYTHIA8_DIRS} PATH_SUFFIXES include include/Pythia8 include/pythia8)

