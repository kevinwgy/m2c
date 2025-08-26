This is a test case that demonstrates the 2D->3D projection feature for fluid-structure coupled simulations.
M2C operates on a 2D-cylindrical mesh, while the structural solver, Aero-S, operates on a 3D mesh.

Command (example): mpiexec -n 31 [path-to-m2c-executible] input.st : -n 1 [path-to-aeros-executible] fem.in 2>&1 | tee log.out

Tested with Aero-S changeset 3152:f484d5c512c8

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

