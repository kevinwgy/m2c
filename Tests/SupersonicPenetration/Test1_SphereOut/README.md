# Overview
An example problem featuring a spherical impactor penetrating a thin-walled rectangular container at supersonic speed. Two different fluid materials inside and outside the container. The embedded boundary method in M2C handles the topological change of the computational domain resulting from contact and penetration.

# Instructions
1. Compile the user-specified impactor dynamics file. (See header of `UserDefinedDynamics.cpp`.

2. Run the simulation: See job submission script `run.st`.

3. Post-process the embedded boundary dynamics for visualization
```
cd results
xp2exo impactor.top impactor.exo impactor_disp.txt
xp2exo target.top target.exo target_disp.txt
```
