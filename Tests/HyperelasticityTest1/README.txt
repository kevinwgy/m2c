-------------
10/31/2023
K.W.
-------------

This test case calls the hyperelasticity models to calculate the principal elastic stresses at 
one or multiple probe locations. There are eight input files in this folder, corresponding to
3D and 2D-Cylindrical domains and in each case, four different hyperelasticity models 
(S.Venant-Kirchhoff, Modified S.Venant-Kirchhoff, Neo-Hookean, and Mooney-Rivlin). 

A simple shear is prescibed. The problem setup and exact solutions are given in KW's notes.

The log files and results (in the "results" folder) are generated using the following version
of the M2C code:
----------------------------------------------------------------------
commit 612a8579f5e525d9f8fa00ca92b3f5cd8719c943
Author: Kevin Guanyuan Wang <kevin.wgy@gmail.com>
Date:   Tue Oct 31 21:30:16 2023 -0400

    Minor bug fixes. Verified probe-output of principal stresses.
----------------------------------------------------------------------
NOTE: Before running the test, the following flag in CMakeLists must be turned on.
--------------------------------------------
add_definitions(-DHYPERELASTICITY_TEST=1)
--------------------------------------------

The run command can be found in the log files. The principal elastic stresses (i.e. outputs)
are written into "probe_stresses_XXX.txt" in the "results" folder.

