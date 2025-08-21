This example demonstrates how to restart a simulation on a different mesh using the 
UserDefinedState feature.

Within this repository, "input.st" is the input file of the original simulation, 
which involves laser-fluid coupling, laser-induced vaporization, and two embedded
boundaries. "IC/solution.vtr" is a solution snapshot obtained from this simulation.
"input_restart.st" is the input file for restarting the simulation. It differs from
the original one in that
(1) The mesh size and resolutions are different.
(2) User-defined initial state is activated through
    UserDefinedStateCalculator = "IC/libUserDefinedState.so";
(3) Laser radiation is turned off. (While laser radiation can be included in the
    restart simulation, the UserDefinedState feature does not currently support
    initializing the laser irradiance field.)

In "IC/UserDefinedState.cpp", function "MyStateCalculator::GetUserDefinedState" has
been implemented. It reads the binary VTR (VTK Rectilinear) solution file ("solution.vtr"),
extracts the original mesh and solution fields, and use trilinear interpolation to project
these solution fields onto the new mesh, which serve as the initial condition of the 
restart simulation.

NOTE: Users should modify this function to fit their specific problem.

To run this example:
- Enter subdirectory "IC".
- Open "CMakeLists.txt", update the path in "target_include_directories". It should point
  to the main repository of the M2C code.
- Compile "UserDefinedState.cpp" as a shared object:
  $cmake .
  $make
  It should generate file "libUserDefinedState.so".
- Return to the simulation folder. Launch the restart simulation
  (e.g., "mpiexec -n [N] [m2c-executible] input_restart.st").
