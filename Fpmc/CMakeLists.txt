include_directories(External ${PROJECT_SOURCE_DIR}/Herwig)
include_directories(interface)

set(fpmc_sources src/fpmc.f src/fpmc_welcome.f src/Fpmc.cc src/FpmcParameters.cc)
if(CERNLIB)
  list(APPEND fpmc_sources src/ffkey.f)
endif()

add_subdirectory(External)

add_library(FPMCDeps SHARED EXCLUDE_FROM_ALL $<TARGET_OBJECTS:CHIDe> $<TARGET_OBJECTS:KMR_lumi> $<TARGET_OBJECTS:KMR_lumi_2> $<TARGET_OBJECTS:excl_aaaa> $<TARGET_OBJECTS:softc> $<TARGET_OBJECTS:miscellaneous>)
target_link_libraries(FPMCDeps comphep)
set_target_properties(FPMCDeps PROPERTIES LINKER_LANGUAGE CXX)

add_library(FPMC SHARED ${fpmc_sources} $<TARGET_OBJECTS:CHIDe> $<TARGET_OBJECTS:KMR_lumi> $<TARGET_OBJECTS:KMR_lumi_2> $<TARGET_OBJECTS:excl_aaaa> $<TARGET_OBJECTS:softc> $<TARGET_OBJECTS:miscellaneous>)
target_link_libraries(FPMC comphep)
set_target_properties(FPMC PROPERTIES LINKER_LANGUAGE CXX)

