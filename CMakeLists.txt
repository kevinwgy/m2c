cmake_minimum_required(VERSION 3.10)

# set the project name
project(M2C VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# use mpicc and mpicxx to compile
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# add path to find_package scripts
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# compiler flags (turn on all the warnings)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-reorder")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# folders for executable and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

#--------------------------------------------------------
# find packages 
if(NOT DEFINED ENV{PETSC_DIR})
  message(FATAL_ERROR "ERROR: Could not find Petsc. Linux environment variable PETSC_DIR needs to be set. Try 'export PETSC_IDR=[path-to-petsc-build-directory]'.")
endif()

find_package(PETSc REQUIRED)
if(NOT "${PETSc_LIBRARY_DIR}" STREQUAL "")
  link_directories("${PETSc_LIBRARY_DIR}")
endif()
if(NOT "${PETSc_INCLUDE_DIR}" STREQUAL "")
  include_directories("${PETSc_INCLUDE_DIR}")
endif()

find_package(Eigen3 3.3 REQUIRED)
if(NOT "${EIGEN3_INCLUDE_DIR}" STREQUAL "")
  message("-- Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
  include_directories("${EIGEN3_INCLUDE_DIR}")
else()
  message(FATAL_ERROR "ERROR: Found Eigen3 but EIGEN3_INCLUDE_DIR undefined. (Try again with -DEIGEN3_INCLUDE_DIR=<provide-path-here>")
endif()

find_package(Boost 1.72 REQUIRED)
if(NOT "${Boost_INCLUDE_DIR}" STREQUAL "")
  include_directories("${Boost_INCLUDE_DIR}")
endif()

#--------------------------------------------------------

#-----------------------------------------
# Special tag(s). Don't turn on unless you know what they mean.
#
#add_definitions(-DLEVELSET_TEST=3)
#add_definitions(-DHYPERELASTICITY_TEST=1)
#-----------------------------------------

# -----------------------------
# for version control
# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
add_custom_target(VersionHeader
  ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/version.cmake
)

ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_SOURCE_DIR}/version.h
           ${CMAKE_SOURCE_DIR}/_version.h
    COMMAND cmake -P
            ${CMAKE_SOURCE_DIR}/version.h)
# -----------------------------

# add include and link directories 
include_directories("${CMAKE_SOURCE_DIR}")
include_directories(${CMAKE_SOURCE_DIR}/GeoTools)
include_directories(${CMAKE_SOURCE_DIR}/MathTools)
link_directories(${CMAKE_SOURCE_DIR}/parser)

# compile the parser
add_custom_target(
   extern_lib
   COMMAND $(MAKE)
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/parser
)

# add the executable
add_executable(m2c 
Main.cpp
IoData.cpp
SpaceVariable.cpp
ConcurrentProgramsHandler.cpp
CommunicationTools.cpp
AerosMessenger.cpp
AerofMessenger.cpp
M2CTwinMessenger.cpp
TriangulatedSurface.cpp
CrackingSurface.cpp
Intersector.cpp
FloodFill.cpp
MultiSurfaceIntersector.cpp
EmbeddedBoundaryOperator.cpp
CustomCommunicator.cpp
ExactRiemannSolverBase.cpp
ExactRiemannSolverInterfaceJump.cpp
MultiPhaseOperator.cpp
SymmetryOperator.cpp
SmoothingOperator.cpp
Interpolator.cpp
GlobalMeshInfo.cpp
NeighborCommunicator.cpp
GradientCalculatorCentral.cpp
GradientCalculatorFD3.cpp
ViscosityOperator.cpp
HeatDiffusionOperator.cpp
ProbeOutput.cpp
EnergyIntegrationOutput.cpp
IntegrationOutput.cpp
PhaseTransitionOutput.cpp
PlaneOutput.cpp
MaterialVolumeOutput.cpp
Output.cpp
SpaceInitializer.cpp
Reconstructor.cpp
SpaceOperator.cpp
IncompressibleOperator.cpp
MeshGenerator.cpp
MeshMatcher.cpp
LevelSetOperator.cpp
LevelSetReinitializer.cpp
LaserAbsorptionSolver.cpp
IonizationOperator.cpp
SahaEquationSolver.cpp
NonIdealSahaEquationSolver.cpp
AtomicIonizationData.cpp
EmbeddedBoundaryFormula.cpp
HyperelasticityFcn.cpp
HyperelasticityFcn2DCyl.cpp
HyperelasticityOperator.cpp
ReferenceMapOperator.cpp
PrescribedMotionOperator.cpp
TimeIntegrator.cpp
TimeIntegratorSemiImp.cpp
SpecialToolsDriver.cpp
DynamicLoadCalculator.cpp
LagrangianOutput.cpp
EOSAnalyzer.cpp
TerminalVisualization.cpp
ClosestTriangle.cpp
SteadyStateOperator.cpp
GhostFluidOperator.cpp
GravityHandler.cpp
LinearOperator.cpp
LinearSystemSolver.cpp
Utils.cpp
MathTools/rbf_interp.cpp
MathTools/polynomial_equations.cpp
MathTools/linear_algebra_2d.cpp
MathTools/linear_algebra_3d.cpp
GeoTools/GeoTools.cpp
GeoTools/BoundingBoxes.cpp
GeoTools/Intersections.cpp)


# link to libraries
target_link_libraries(m2c petsc mpi parser)
target_link_libraries(m2c ${CMAKE_DL_LIBS}) #linking to the dl library (-ldl)
add_dependencies(m2c extern_lib)
add_dependencies(m2c VersionHeader)
