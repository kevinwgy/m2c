This is a fluid-structure coupled simulation of near-field underwater explosion (UNDEX).
It is similar to the problem presented in Wentao Ma et al., IJSS (2022), except that
the fluid and structural meshes are both coarser here. 

tinkercliffs_sbatch.sh --- job submission script 
input.st --- M2C input file
fem.in, mesh.include --- Aero-S input file


---------------------------------------
Version of Aero-S used in this example
(tip of hg log)
---------------------------------------

changeset:   3152:f484d5c512c8
tag:         tip
user:        Kevin G. Wang <kevinwgy@vt.edu>
date:        Tue Aug 27 09:20:59 2019 -0400
files:       Dist.d/DistDom.C Driver.d/BinaryOutputInstance.C Driver.d/DecDomainImpl.h Driver.d/NLStatic.C Driver.d/OpMake.C Driver.d/SOpt.C
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
Implemented output of a residual-based, a posteriori, relative error indicator for assessing the accuracy of a PROM; activate using keyword errindic under OUTPUT for dynamic analysis with srom/READMODE.


changeset:   3149:aa79fe083700
user:        Philip Avery <pavery@stanford.edu>
date:        Mon Jul 29 22:59:24 2019 -0700
files:       Driver.d/Header.h
description:
Fixed bug in output file header.
