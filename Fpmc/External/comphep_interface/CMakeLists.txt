set(processes sqme_aaaa sqme_aazz sqme_aaww)

#----- build all dependencies

include_directories(commons)

set(comphep_sources "")
foreach(_proc ${processes})
  include_directories(${_proc}/comphep/include)
  file(GLOB sources ${_proc}/src/*.cpp)
  list(APPEND comphep_sources ${sources})
endforeach()

#----- build a big library out of all sub-processes

add_library(comphep OBJECT ${comphep_sources})
set_target_properties(comphep PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

