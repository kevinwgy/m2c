cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(UserDefinedState)

find_package(VTK COMPONENTS 
  CommonCore
  IOXML
  CommonDataModel   
  FiltersSources
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "UserDefinedState: Unable to find the VTK build folder.")
endif()

add_library(UserDefinedState SHARED UserDefinedState.cpp)

# ----------------------------------------------------------------
# [ATTENTION] Provide the path to the main folder of the M2C code
target_include_directories(UserDefinedState PUBLIC ../../..)
# ----------------------------------------------------------------

target_link_libraries(UserDefinedState PUBLIC ${VTK_LIBRARIES})
# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS UserDefinedState
  MODULES ${VTK_LIBRARIES}
)
