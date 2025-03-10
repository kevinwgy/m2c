This folder contains three test cases for the level set method. The fluid governing equations are NOT solved. Instead, the velocity field is pre-specified.

NOTE: To run these cases, the environment variable "LEVELSET_TEST" (in CMakeLists.txt) must be set:
add_definitions(-DLEVELSET_TEST=1) ~~~ SlottedDisk
add_definitions(-DLEVELSET_TEST=2) ~~~ VortexDeformation
add_definitions(-DLEVELSET_TEST=3) ~~~ MergingCircles
