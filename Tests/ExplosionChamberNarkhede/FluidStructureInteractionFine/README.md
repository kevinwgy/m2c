# Overview
This is an example case that simulates the dynamics of a containment
structure that is subjected to an internal explosion. The initial
conditions within the confinement structure are provided through
the text file `SphericalShock.txt`. 

For details regarding the simulation setup, refer
`Fluid–structure coupled simulation framework for lightweight explosion
containment structures under large deformations` by Narkhede et. al.


# Structural dynamics solver
This test case utilizes the Aero-S solver as the structural dynamics solver. It has been tested with Aero-S changeset 3152:f484d5c512c8

changeset:   3152:f484d5c512c8
tag:         tip
user:        Kevin G. Wang <kevinwgy@vt.edu>
date:        Tue Aug 27 09:20:59 2019 -0400
files:       Dist.d/DistDom.C Driver.d/BinaryOutputInstance.C Driver.d/DecDomainImpl.h Driver.d/NLStatic.C Driver.d/OpM
ake.C Driver.d/SOpt.C
description:
Enabled output (binary & ascii) for principal axes for certain types of simulations.

changeset:   3151:aeb023d3c090
user:        Kevin G. Wang <kevinwgy@vt.edu>
date:        Sun Aug 25 07:42:51 2019 -0400
files:       Element.d/BelytschkoTsayShell.d/BelytschkoTsayShell.C
description:
Clarified some ambiguities about Vector in BelytschkoTsayShell.C (now using ::Vector).

changeset:   3150:4823b184e11d
user:        Philip Avery <pavery@stanford.edu>
date:        Tue Jul 30 09:08:40 2019 -0700
files:       Driver.d/Dynam.C Driver.d/Header.h Parser.d/lexer.l Problems.d/ModalBase.C Utils.d/OutputInfo.h
description:
Implemented output of a residual-based, a posteriori, relative error indicator for assessing the accuracy of a PROM; ac
tivate using keyword errindic under OUTPUT for dynamic analysis with srom/READMODE.

# Instructions on post-processing Aero-S results
In this case, structural results are written in binary files (in `fem.st`: `BINARYOUTPUT On`). It saves computation time, especially when the structural model is large. Post-processing the results for visualization takes a few steps:

1. Create a structural mesh file with only `NODES` and `TOPOLOGY`, with an `END` at the end of the file. See, e.g., `mesh.include` in sub-directory `aeros_postpro`.

2. Create binary input files using the `Sower` tool (but only used for postprocessing results). For example,
```
cd aeros_postpro
sower -struct -mesh mesh.include -dec ../ExampleCapsule.optDec -cpu 16 -cluster 1
```
Note that the mesh partition file `ExampleCapsule.optDec` will be automatically generated after launching the simulation.
Now, there should be `OUTPUT.16cpu`, `OUTPUT.con`, `OUTPUT.dec1`, and `OUTPOUT.msh1` in the sub-directory.

3. Convert binary results into ASCII files, again using `Sower`. For example,
```
cd aeros_postpro
sower -struct -merge -con OUTPUT.con -mesh OUTPUT.msh -result ../results/disp 
sower -struct -merge -con OUTPUT.con -mesh OUTPUT.msh -result ../results/epstrain
```
These commands should generate files `disp.xpost` and `epstrain.xpost` in sub-directory `results` (instead of the local directory).

4. (Optional) Combine structural mesh and results (both in ASCII format) into an `.exo` file (Exodus II format) for visualization, using the `xp2exo` tool.
```
[Aero-S executable] -t fem.in
cd results
xp2exo ../ExampleCapsule.top capsule.exo disp.xpost epstrain.xpost
```
The first command takes in the input file `fem.in`, and outputs the mesh file `ExampleCapsule.top`. Running the `xp2exo` command should create ``capsule.exo``, which can be visualized by Paraview and other visualization software.
