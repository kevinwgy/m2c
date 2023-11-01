/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <Utils.h>
#include <IoData.h>
#include <parser/Assigner.h>
#include <parser/Dictionary.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cfloat>
#include <climits>
#include <cmath>
#include <unistd.h>
#include <bits/stdc++.h> //INT_MAX
//#include <dlfcn.h>
using namespace std;

double avogadro_number = 6.02214076e23;

RootClassAssigner *nullAssigner = new RootClassAssigner;

//------------------------------------------------------------------------------

StateVariable::StateVariable()
{

  materialid  = 0;
  density     = 1.0e-6;
  velocity_x  = 0.0;
  velocity_y  = 0.0;
  velocity_z  = 0.0;
  pressure    = 0.0;
  temperature = 0.0;
  internal_energy_per_mass = 0.0;

}

//------------------------------------------------------------------------------

void StateVariable::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 8, father);

  new ClassInt<StateVariable>(ca, "MaterialID", this, &StateVariable::materialid);
  new ClassDouble<StateVariable>(ca, "Density", this, &StateVariable::density);
  new ClassDouble<StateVariable>(ca, "VelocityX", this, &StateVariable::velocity_x);
  new ClassDouble<StateVariable>(ca, "VelocityY", this, &StateVariable::velocity_y);
  new ClassDouble<StateVariable>(ca, "VelocityZ", this, &StateVariable::velocity_z);
  new ClassDouble<StateVariable>(ca, "Pressure", this, &StateVariable::pressure);
  new ClassDouble<StateVariable>(ca, "Temperature", this, &StateVariable::temperature);
  new ClassDouble<StateVariable>(ca, "InternalEnergyPerUnitMass", this, &StateVariable::internal_energy_per_mass);

}

//------------------------------------------------------------------------------

PointData::PointData()
{
  x  = 0.0;
  y  = 0.0;
  z  = 0.0;

  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *PointData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 6, nullAssigner);

  new ClassDouble<PointData>
    (ca, "X", this, &PointData::x);
  new ClassDouble<PointData>
    (ca, "Y", this, &PointData::y);
  new ClassDouble<PointData>
    (ca, "Z", this, &PointData::z);

  initialConditions.setup("InitialState", ca);

  new ClassToken<PointData> (ca, "Inclusion", this,
     reinterpret_cast<int PointData::*>(&PointData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<PointData>(ca, "OperationOrder", this, &PointData::order);

  return ca;
}

//------------------------------------------------------------------------------


PlaneData::PlaneData()
{

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  nx     = 0.0;
  ny     = 0.0;
  nz     = 0.0;

  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *PlaneData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 9, nullAssigner);

  new ClassDouble<PlaneData> (ca, "Point_x", this, &PlaneData::cen_x);
  new ClassDouble<PlaneData> (ca, "Point_y", this, &PlaneData::cen_y);
  new ClassDouble<PlaneData> (ca, "Point_z", this, &PlaneData::cen_z);
  new ClassDouble<PlaneData> (ca, "Normal_x", this, &PlaneData::nx);
  new ClassDouble<PlaneData> (ca, "Normal_y", this, &PlaneData::ny);
  new ClassDouble<PlaneData> (ca, "Normal_z", this, &PlaneData::nz);

  initialConditions.setup("InitialState", ca);

  new ClassToken<PlaneData> (ca, "Inclusion", this,
     reinterpret_cast<int PlaneData::*>(&PlaneData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<PlaneData>(ca, "OperationOrder", this, &PlaneData::order);

  return ca;
}

//------------------------------------------------------------------------------

SphereData::SphereData()
{
  
  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  radius = -1.0;

  side = INTERIOR;
  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *SphereData::getAssigner()
{
  
  ClassAssigner *ca = new ClassAssigner("normal", 8, nullAssigner);
  
  new ClassDouble<SphereData> (ca, "Center_x", this, &SphereData::cen_x);
  new ClassDouble<SphereData> (ca, "Center_y", this, &SphereData::cen_y);
  new ClassDouble<SphereData> (ca, "Center_z", this, &SphereData::cen_z);
  new ClassDouble<SphereData> (ca, "Radius", this, &SphereData::radius);
  
  new ClassToken<SphereData> (ca, "Side", this,
     reinterpret_cast<int SphereData::*>(&SphereData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<SphereData> (ca, "Inclusion", this,
     reinterpret_cast<int SphereData::*>(&SphereData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<SphereData>(ca, "OperationOrder", this, &SphereData::order);

  initialConditions.setup("InitialState", ca);
  
  return ca;
}

//------------------------------------------------------------------------------

ParallelepipedData::ParallelepipedData()
{

  x0 = y0 = z0 = 0.0;
  ax = ay = az = 0.0;
  bx = by = bz = 0.0;
  cx = cy = cz = 0.0;

  side = INTERIOR;
  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *ParallelepipedData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 16, nullAssigner);

  new ClassDouble<ParallelepipedData> (ca, "Px", this, &ParallelepipedData::x0);
  new ClassDouble<ParallelepipedData> (ca, "Py", this, &ParallelepipedData::y0);
  new ClassDouble<ParallelepipedData> (ca, "Pz", this, &ParallelepipedData::z0);

  new ClassDouble<ParallelepipedData> (ca, "Ax", this, &ParallelepipedData::ax);
  new ClassDouble<ParallelepipedData> (ca, "Ay", this, &ParallelepipedData::ay);
  new ClassDouble<ParallelepipedData> (ca, "Az", this, &ParallelepipedData::az);

  new ClassDouble<ParallelepipedData> (ca, "Bx", this, &ParallelepipedData::bx);
  new ClassDouble<ParallelepipedData> (ca, "By", this, &ParallelepipedData::by);
  new ClassDouble<ParallelepipedData> (ca, "Bz", this, &ParallelepipedData::bz);

  new ClassDouble<ParallelepipedData> (ca, "Cx", this, &ParallelepipedData::cx);
  new ClassDouble<ParallelepipedData> (ca, "Cy", this, &ParallelepipedData::cy);
  new ClassDouble<ParallelepipedData> (ca, "Cz", this, &ParallelepipedData::cz);

  new ClassToken<ParallelepipedData> (ca, "Side", this,
     reinterpret_cast<int ParallelepipedData::*>(&ParallelepipedData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<ParallelepipedData> (ca, "Inclusion", this,
     reinterpret_cast<int ParallelepipedData::*>(&ParallelepipedData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<ParallelepipedData>(ca, "OperationOrder", this, &ParallelepipedData::order);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

SpheroidData::SpheroidData()
{

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;

  axis_x = 0.0;
  axis_y = 0.0;
  axis_z = 0.0;

  semi_length = 0.0;
  radius = 0.0;

  side = INTERIOR;
  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *SpheroidData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 12, nullAssigner);

  new ClassDouble<SpheroidData> (ca, "Center_x", this, &SpheroidData::cen_x);
  new ClassDouble<SpheroidData> (ca, "Center_y", this, &SpheroidData::cen_y);
  new ClassDouble<SpheroidData> (ca, "Center_z", this, &SpheroidData::cen_z);
  new ClassDouble<SpheroidData> (ca, "Axis_x", this, &SpheroidData::axis_x);
  new ClassDouble<SpheroidData> (ca, "Axis_y", this, &SpheroidData::axis_y);
  new ClassDouble<SpheroidData> (ca, "Axis_z", this, &SpheroidData::axis_z);
  new ClassDouble<SpheroidData> (ca, "SemiLength", this, &SpheroidData::semi_length);
  new ClassDouble<SpheroidData> (ca, "Radius", this, &SpheroidData::radius);

  new ClassToken<SpheroidData> (ca, "Side", this,
     reinterpret_cast<int SpheroidData::*>(&SpheroidData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<SpheroidData> (ca, "Inclusion", this,
     reinterpret_cast<int SpheroidData::*>(&SpheroidData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<SpheroidData>(ca, "OperationOrder", this, &SpheroidData::order);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

CylinderConeData::CylinderConeData() {

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  nx     = 0.0;
  ny     = 0.0;
  nz     = 1.0;

  r = 1.0;
  L = 0.0;

  cone_height = 0.0;
  opening_angle_degrees = 45.0;

  side = INTERIOR;
  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *CylinderConeData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 14, nullAssigner);

  new ClassDouble<CylinderConeData> (ca, "Axis_x", this, &CylinderConeData::nx);
  new ClassDouble<CylinderConeData> (ca, "Axis_y", this, &CylinderConeData::ny);
  new ClassDouble<CylinderConeData> (ca, "Axis_z", this, &CylinderConeData::nz);
  new ClassDouble<CylinderConeData> (ca, "BaseCenter_x", this, &CylinderConeData::cen_x);
  new ClassDouble<CylinderConeData> (ca, "BaseCenter_y", this, &CylinderConeData::cen_y);
  new ClassDouble<CylinderConeData> (ca, "BaseCenter_z", this, &CylinderConeData::cen_z);
  new ClassDouble<CylinderConeData> (ca, "CylinderRadius", this, &CylinderConeData::r);
  new ClassDouble<CylinderConeData> (ca, "CylinderHeight", this, &CylinderConeData::L);

  new ClassDouble<CylinderConeData> (ca, "ConeOpeningAngleInDegrees", this, &CylinderConeData::opening_angle_degrees);
  new ClassDouble<CylinderConeData> (ca, "ConeHeight", this, &CylinderConeData::cone_height);

  new ClassToken<CylinderConeData> (ca, "Side", this,
     reinterpret_cast<int CylinderConeData::*>(&CylinderConeData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<CylinderConeData> (ca, "Inclusion", this,
     reinterpret_cast<int CylinderConeData::*>(&CylinderConeData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<CylinderConeData>(ca, "OperationOrder", this, &CylinderConeData::order);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

CylinderSphereData::CylinderSphereData() {

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  nx     = 0.0;
  ny     = 0.0;
  nz     = 1.0;

  r = 1.0;
  L = 0.0;

  front_cap = Off;
  back_cap = Off;

  side = INTERIOR;
  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *CylinderSphereData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 14, nullAssigner);

  new ClassDouble<CylinderSphereData> (ca, "Axis_x", this, &CylinderSphereData::nx);
  new ClassDouble<CylinderSphereData> (ca, "Axis_y", this, &CylinderSphereData::ny);
  new ClassDouble<CylinderSphereData> (ca, "Axis_z", this, &CylinderSphereData::nz);
  new ClassDouble<CylinderSphereData> (ca, "BaseCenter_x", this, &CylinderSphereData::cen_x);
  new ClassDouble<CylinderSphereData> (ca, "BaseCenter_y", this, &CylinderSphereData::cen_y);
  new ClassDouble<CylinderSphereData> (ca, "BaseCenter_z", this, &CylinderSphereData::cen_z);
  new ClassDouble<CylinderSphereData> (ca, "CylinderRadius", this, &CylinderSphereData::r);
  new ClassDouble<CylinderSphereData> (ca, "CylinderHeight", this, &CylinderSphereData::L);
  new ClassToken<CylinderSphereData> (ca, "FrontSphericalCap", this,
     reinterpret_cast<int CylinderSphereData::*>(&CylinderSphereData::front_cap), 2,
     "Off", 0, "On", 1);
  new ClassToken<CylinderSphereData> (ca, "BackSphericalCap", this,
     reinterpret_cast<int CylinderSphereData::*>(&CylinderSphereData::back_cap), 2,
     "Off", 0, "On", 1);

  new ClassToken<CylinderSphereData> (ca, "Side", this,
     reinterpret_cast<int CylinderSphereData::*>(&CylinderSphereData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<CylinderSphereData> (ca, "Inclusion", this,
     reinterpret_cast<int CylinderSphereData::*>(&CylinderSphereData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<CylinderSphereData>(ca, "OperationOrder", this, &CylinderSphereData::order);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

UserSpecifiedEnclosureData::UserSpecifiedEnclosureData()
{
  surface_filename = "";    
  surface_thickness = 1.0e-8;

  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *UserSpecifiedEnclosureData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 4, nullAssigner);

  new ClassStr<UserSpecifiedEnclosureData>(ca, "SurfaceMeshFile", this,
          &UserSpecifiedEnclosureData::surface_filename);

  new ClassDouble<UserSpecifiedEnclosureData>(ca, "SurfaceThickness", this, 
          &UserSpecifiedEnclosureData::surface_thickness);

  initialConditions.setup("InitialState", ca);

  new ClassToken<UserSpecifiedEnclosureData> (ca, "Inclusion", this,
     reinterpret_cast<int UserSpecifiedEnclosureData::*>(&UserSpecifiedEnclosureData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<UserSpecifiedEnclosureData>(ca, "OperationOrder", this, &UserSpecifiedEnclosureData::order);

  return ca;
}

//------------------------------------------------------------------------------

void MultiInitialConditionsData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 8, father);
  pointMap.setup("Point", ca);
  planeMap.setup("Plane", ca);
  sphereMap.setup("Sphere", ca);
  parallelepipedMap.setup("Parallelepiped", ca);
  spheroidMap.setup("Spheroid", ca);
  cylinderconeMap.setup("CylinderAndCone", ca);
  cylindersphereMap.setup("CylinderWithSphericalCaps", ca);
  enclosureMap.setup("ArbitraryEnclosure", ca);
}

//------------------------------------------------------------------------------

MeshResolution1DPointData::MeshResolution1DPointData()
{
  coord = 0.0;
  h     = 0.0;
}

//------------------------------------------------------------------------------

Assigner *MeshResolution1DPointData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 2, nullAssigner);

  new ClassDouble<MeshResolution1DPointData>
    (ca, "Coordinate", this, &MeshResolution1DPointData::coord);
  new ClassDouble<MeshResolution1DPointData>
    (ca, "CellWidth", this, &MeshResolution1DPointData::h);

  return ca;
}

//------------------------------------------------------------------------------

MeshData::MeshData()
{
  type = THREEDIMENSIONAL;
  x0 = 0.0;
  xmax = 1.0;
  y0 = 0.0;
  ymax = 1.0;
  z0 = 0.0;
  zmax = 1.0;
  Nx = -1;
  Ny = -1;
  Nz = -1;

  bc_x0   = NONE;
  bc_xmax = NONE;
  bc_y0   = NONE;
  bc_ymax = NONE;
  bc_z0   = NONE;
  bc_zmax = NONE;
}

//------------------------------------------------------------------------------

void MeshData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 20, father);

  new ClassToken<MeshData>(ca, "Type", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::type), 3,
                               "ThreeDimensional", 0, "Spherical", 1, "Cylindrical", 2);
  new ClassDouble<MeshData>(ca, "X0", this, &MeshData::x0);
  new ClassDouble<MeshData>(ca, "Xmax", this, &MeshData::xmax);
  new ClassDouble<MeshData>(ca, "Y0", this, &MeshData::y0);
  new ClassDouble<MeshData>(ca, "Ymax", this, &MeshData::ymax);
  new ClassDouble<MeshData>(ca, "Z0", this, &MeshData::z0);
  new ClassDouble<MeshData>(ca, "Zmax", this, &MeshData::zmax);
  new ClassInt<MeshData>(ca, "NumberOfCellsX", this, &MeshData::Nx);
  new ClassInt<MeshData>(ca, "NumberOfCellsY", this, &MeshData::Ny);
  new ClassInt<MeshData>(ca, "NumberOfCellsZ", this, &MeshData::Nz);

  xpoints_map.setup("ControlPointX", ca);
  ypoints_map.setup("ControlPointY", ca);
  zpoints_map.setup("ControlPointZ", ca);

  // Inside the code: Farfield0 = Farfield = Inlet, Farfield1 = Outlet
  new ClassToken<MeshData>(ca, "BoundaryConditionX0", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_x0), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
  new ClassToken<MeshData>(ca, "BoundaryConditionXmax", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_xmax), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
  new ClassToken<MeshData>(ca, "BoundaryConditionY0", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_y0), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
  new ClassToken<MeshData>(ca, "BoundaryConditionYmax", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_ymax), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
  new ClassToken<MeshData>(ca, "BoundaryConditionZ0", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_z0), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
  new ClassToken<MeshData>(ca, "BoundaryConditionZmax", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_zmax), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
 } 

//------------------------------------------------------------------------------

void MeshData::check()
{
  if(type == SPHERICAL) {
    if(Ny != 1 || Nz != 1) {
      print_error("*** Error: For a domain w/ spherical symmetry, the mesh should "
                  "have 1 cell in y- and z-dirs (%d, %d).\n", Ny, Nz); 
      exit_mpi();
    }
    if(x0<0.0) {
      print_error("*** Error: For a domain w/ spherical symmetry, x0 must be nonnegative (%e).\n",
                  x0);
      exit_mpi();
    }
  }
  else if(type == CYLINDRICAL) {
    if(Nz != 1) {
      print_error("*** Error: For a domain w/ cylindrical symmetry, the mesh should "
                  "have 1 cell in z-dir (%d).\n", Nz); 
      exit_mpi();
    }
    if(y0<0.0) {
      print_error("*** Error: For a domain w/ cylindrical symmetry, y0 must be nonnegative (%e).\n",
                  y0);
      exit_mpi();
    }
  }
}

//------------------------------------------------------------------------------

StiffenedGasModelData::StiffenedGasModelData()
{

  specificHeatRatio = 1.4;
  pressureConstant = 0.0;
  enthalpyConstant = 0.0;

  cv = 0.0;
  T0 = 0.0;
  e0 = 0.0;

  cp = 0.0;
  h0 = 0.0;

  rho0 = 0.0;

}

//------------------------------------------------------------------------------

void StiffenedGasModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 9, father);

  new ClassDouble<StiffenedGasModelData>(ca, "SpecificHeatRatio", this,
                                &StiffenedGasModelData::specificHeatRatio);
  new ClassDouble<StiffenedGasModelData>(ca, "PressureConstant", this,
                                &StiffenedGasModelData::pressureConstant);
  new ClassDouble<StiffenedGasModelData>(ca, "EnthalpyConstant", this,
                                &StiffenedGasModelData::enthalpyConstant);

  new ClassDouble<StiffenedGasModelData>(ca, "SpecificHeatAtConstantVolume", this,
                                &StiffenedGasModelData::cv);
  new ClassDouble<StiffenedGasModelData>(ca, "ReferenceTemperature", this,
                                &StiffenedGasModelData::T0);
  new ClassDouble<StiffenedGasModelData>(ca, "ReferenceSpecificInternalEnergy", this,
                                &StiffenedGasModelData::e0);

  new ClassDouble<StiffenedGasModelData>(ca, "SpecificHeatAtConstantPressure", this,
                                &StiffenedGasModelData::cp);
  new ClassDouble<StiffenedGasModelData>(ca, "ReferenceSpecificEnthalpy", this,
                                &StiffenedGasModelData::h0);

  new ClassDouble<StiffenedGasModelData>(ca, "ReferenceDensity", this,
                                &StiffenedGasModelData::rho0);

}

//------------------------------------------------------------------------------

NobleAbelStiffenedGasModelData::NobleAbelStiffenedGasModelData()
{

  specificHeatRatio = 1.4;
  pressureConstant  = 0.0;
  volumeConstant    = 0.0;
  energyConstant    = 0.0;
  entropyConstant   = 0.0;

  cv = 0.0;

}

//------------------------------------------------------------------------------

void NobleAbelStiffenedGasModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassDouble<NobleAbelStiffenedGasModelData>(ca, "SpecificHeatRatio", this,
                      &NobleAbelStiffenedGasModelData::specificHeatRatio);
  new ClassDouble<NobleAbelStiffenedGasModelData>(ca, "PressureConstant", this,
                      &NobleAbelStiffenedGasModelData::pressureConstant);
  new ClassDouble<NobleAbelStiffenedGasModelData>(ca, "VolumeConstant", this,
                      &NobleAbelStiffenedGasModelData::volumeConstant);
  new ClassDouble<NobleAbelStiffenedGasModelData>(ca, "EnergyConstant", this,
                      &NobleAbelStiffenedGasModelData::energyConstant);
  new ClassDouble<NobleAbelStiffenedGasModelData>(ca, "EntropyConstant", this,
                      &NobleAbelStiffenedGasModelData::entropyConstant);

  new ClassDouble<NobleAbelStiffenedGasModelData>(ca, "SpecificHeatAtConstantVolume", this,
                      &NobleAbelStiffenedGasModelData::cv);

}


//------------------------------------------------------------------------------

MieGruneisenModelData::MieGruneisenModelData()
{
  // default values are for copper
  // These values are taken from Wikipedia
  // (https://en.wikipedia.org/wiki/Mie%E2%80%93Gr%C3%BCneisen_equation_of_state),
  // which cites two papers: Mitchell and Nellis (1981) "Shock compression of
  // aluminum, copper, and tantalum", Journal of Applied Physics, and MacDonald and
  // MacDonald (1981) "Thermodynamic properties of fcc metals at high temperaures",
  // Physical Review B.

  rho0 = 8.96e-3;       // unit: g/mm3
  c0 = 3.933e6;         // unit: mm/s
  Gamma0 = 1.99;        // non-dimensional
  s = 1.5;              // non-dimensional
  e0 = 0.0;

  cv = 0.0; //3.90e8;   // unit: mm2/(s2.K)
  cp = 0.0;
  h0 = 0.0;
  T0 = 0.0;
}

//------------------------------------------------------------------------------

void MieGruneisenModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 9, father);

  new ClassDouble<MieGruneisenModelData>(ca, "ReferenceDensity", this,
                                         &MieGruneisenModelData::rho0);
  new ClassDouble<MieGruneisenModelData>(ca, "SpecificHeatAtConstantVolume", this,
                                         &MieGruneisenModelData::cv);
  new ClassDouble<MieGruneisenModelData>(ca, "BulkSpeedOfSound", this,
                                         &MieGruneisenModelData::c0);
  new ClassDouble<MieGruneisenModelData>(ca, "HugoniotSlope", this,
                                         &MieGruneisenModelData::s);
  new ClassDouble<MieGruneisenModelData>(ca, "ReferenceGamma", this,
                                         &MieGruneisenModelData::Gamma0);
  new ClassDouble<MieGruneisenModelData>(ca, "ReferenceSpecificInternalEnergy", this,
                                         &MieGruneisenModelData::e0);
  new ClassDouble<MieGruneisenModelData>(ca, "SpecificHeatAtConstantPressure", this,
                                         &MieGruneisenModelData::cp);
  new ClassDouble<MieGruneisenModelData>(ca, "ReferenceSpecificEnthalpy", this,
                                         &MieGruneisenModelData::h0);
  new ClassDouble<MieGruneisenModelData>(ca, "ReferenceTemperature", this,
                                         &MieGruneisenModelData::T0);

}

//------------------------------------------------------------------------------

ExtendedMieGruneisenModelData::ExtendedMieGruneisenModelData() 
{
  // default values are for copper
  // These values are taken from Wikipedia 
  // (https://en.wikipedia.org/wiki/Mie%E2%80%93Gr%C3%BCneisen_equation_of_state),
  // which cites two papers: Mitchell and Nellis (1981) "Shock compression of
  // aluminum, copper, and tantalum", Journal of Applied Physics, and MacDonald and
  // MacDonald (1981) "Thermodynamic properties of fcc metals at high temperaures",
  // Physical Review B.

  rho0 = 8.96e-3;       // unit: g/mm3
  c0 = 3.933e6;         // unit: mm/s
  Gamma0 = 1.99;        // non-dimensional
  s = 1.5;              // non-dimensional
  e0 = 0.0;         

  eta_min = -DBL_MAX;   // non-dimensional (volumetric strain)

  Tlaw = ORIGINAL_CV; 

  cv = 0.0; //3.90e8;   // unit: mm2/(s2.K)
  cp = 0.0;
  h0 = 0.0;
  T0 = 0.0;
}

//------------------------------------------------------------------------------

void ExtendedMieGruneisenModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 11, father);

  new ClassDouble<ExtendedMieGruneisenModelData>(ca, "ReferenceDensity", this, 
                                         &ExtendedMieGruneisenModelData::rho0);
  new ClassDouble<ExtendedMieGruneisenModelData>(ca, "SpecificHeatAtConstantVolume", this, 
                                         &ExtendedMieGruneisenModelData::cv);
  new ClassDouble<ExtendedMieGruneisenModelData>(ca, "BulkSpeedOfSound", this, 
                                         &ExtendedMieGruneisenModelData::c0);
  new ClassDouble<ExtendedMieGruneisenModelData>(ca, "HugoniotSlope", this, 
                                         &ExtendedMieGruneisenModelData::s);
  new ClassDouble<ExtendedMieGruneisenModelData>(ca, "ReferenceGamma", this, 
                                         &ExtendedMieGruneisenModelData::Gamma0);
  new ClassDouble<ExtendedMieGruneisenModelData>(ca, "ReferenceSpecificInternalEnergy", this, 
                                         &ExtendedMieGruneisenModelData::e0);
  new ClassDouble<ExtendedMieGruneisenModelData>(ca, "SpecificHeatAtConstantPressure", this,
                                         &ExtendedMieGruneisenModelData::cp);
  new ClassDouble<ExtendedMieGruneisenModelData>(ca, "ReferenceSpecificEnthalpy", this,
                                         &ExtendedMieGruneisenModelData::h0);
  new ClassDouble<ExtendedMieGruneisenModelData>(ca, "ReferenceTemperature", this,
                                         &ExtendedMieGruneisenModelData::T0);

  new ClassToken<ExtendedMieGruneisenModelData>(ca, "TemperatureLaw", this,
                 reinterpret_cast<int ExtendedMieGruneisenModelData::*>
                 (&ExtendedMieGruneisenModelData::Tlaw), 3,
                 "OriginalCv", 0, "SimplifiedCv", 1, "SimplifiedCp", 2);

  new ClassDouble<ExtendedMieGruneisenModelData>(ca, "VolumetricStrainBreak", this,
                                         &ExtendedMieGruneisenModelData::eta_min);
  
}

//------------------------------------------------------------------------------

TillotsonModelData::TillotsonModelData()
{
  // default values are for water. See Aaron Brundage (2013), Table 1

  rho0 = 0.998e-3;     //unit: g/mm^3
  e0   = 7.0e12;       //unit: mm^2/s^2 = 1.0e-9 J/g = 1.0e-2 erg/g
  a = 0.7;
  b = 0.15;
  A = 2.18e9;          //unit: Pa
  B = 1.325e10;        //unit: Pa
  alpha = 10;
  beta  = 5;

  rhoIV = 0.958e-3;    //unit: g/mm^3
  eIV   = 4.19e11;     //unit: mm^2/s^2
  eCV   = 2.5e12;      //unit: mm^2/s^2

  cv    = 0.0; //3.69e9;      //unit: mm^2/(s^2.K)
  T0    = 0.0;         //unit: K
  temperature_depends_on_density = NO;

  cp    = 0.0;
  h0    = 0.0;
}

//------------------------------------------------------------------------------

void TillotsonModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 16, father);

  new ClassDouble<TillotsonModelData>(ca, "ReferenceDensity", this, &TillotsonModelData::rho0);
  new ClassDouble<TillotsonModelData>(ca, "ReferenceSpecificInternalEnergy", this, &TillotsonModelData::e0);
  new ClassDouble<TillotsonModelData>(ca, "a", this, &TillotsonModelData::a);
  new ClassDouble<TillotsonModelData>(ca, "b", this, &TillotsonModelData::b);
  new ClassDouble<TillotsonModelData>(ca, "A", this, &TillotsonModelData::A);
  new ClassDouble<TillotsonModelData>(ca, "B", this, &TillotsonModelData::B);
  new ClassDouble<TillotsonModelData>(ca, "Alpha", this, &TillotsonModelData::alpha);
  new ClassDouble<TillotsonModelData>(ca, "Beta", this, &TillotsonModelData::beta);
 
  new ClassDouble<TillotsonModelData>(ca, "IncipientVaporizationDensity", this, &TillotsonModelData::rhoIV);
  new ClassDouble<TillotsonModelData>(ca, "IncipientVaporizationSpecificInternalEnergy", this, &TillotsonModelData::eIV);
  new ClassDouble<TillotsonModelData>(ca, "CompleteVaporizationSpecificInternalEnergy", this, &TillotsonModelData::eCV);
 
  

  new ClassDouble<TillotsonModelData>(ca, "SpecificHeatAtConstantVolume", this, &TillotsonModelData::cv);
  new ClassDouble<TillotsonModelData>(ca, "ReferenceTemperature", this, &TillotsonModelData::T0);

  new ClassToken<TillotsonModelData>(ca, "TemperatureDependsOnDensity", this,
                 reinterpret_cast<int TillotsonModelData::*>(&TillotsonModelData::temperature_depends_on_density), 2,
                 "No", 0, "Yes", 1);

  new ClassDouble<TillotsonModelData>(ca, "SpecificHeatAtConstantPressure", this,
                                      &TillotsonModelData::cp);
  new ClassDouble<TillotsonModelData>(ca, "ReferenceSpecificEnthalpy", this,
                                      &TillotsonModelData::h0);

}

//------------------------------------------------------------------------------

JonesWilkinsLeeModelData::JonesWilkinsLeeModelData() 
{
  //default values are for the gaseous products of TNT (Page 21 of Arthur Rallu's thesis)
  omega = 0.28;     //nondimensional
  A1    = 3.712e11; //Pa
  A2    = 3.23e9;   //Pa
  R1    = 4.15;     //nondimensional
  R2    = 0.95;     //nondimensional
  rho0  = 1.63e-3;  //g/mm3

}

//------------------------------------------------------------------------------

void JonesWilkinsLeeModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassDouble<JonesWilkinsLeeModelData>(ca, "Omega", this, 
                                            &JonesWilkinsLeeModelData::omega);
  new ClassDouble<JonesWilkinsLeeModelData>(ca, "A1", this, 
                                            &JonesWilkinsLeeModelData::A1);
  new ClassDouble<JonesWilkinsLeeModelData>(ca, "A2", this, 
                                            &JonesWilkinsLeeModelData::A2);
  new ClassDouble<JonesWilkinsLeeModelData>(ca, "R1", this, 
                                            &JonesWilkinsLeeModelData::R1);
  new ClassDouble<JonesWilkinsLeeModelData>(ca, "R2", this, 
                                            &JonesWilkinsLeeModelData::R2);
  new ClassDouble<JonesWilkinsLeeModelData>(ca, "Rho0", this, 
                                            &JonesWilkinsLeeModelData::rho0);

}

//------------------------------------------------------------------------------

ANEOSBirchMurnaghanDebyeModelData::ANEOSBirchMurnaghanDebyeModelData() 
{
  //Default values are for Copper (J.J. Sanchez, CMAME 2021)
  zeroKelvinDensity = 0.00909; //g/mm3
  b0 = 138.6e9; //Pa
  b0prime = 5.24; //non-D
  delta_e = 0.0;
  molar_mass = 63.546; //g/mol
  T0 = 343.0; //Kelvin
  e0 = 0.0;
  Gamma0 = 1.975; //non-D
  rho0 = 0.00896; //g/mm3

  boltzmann_constant = 1.38064852e-14; //unit: (mm^2).g/(s^2*K)  (dim: [energy]/[temperature])

  debye_evaluation = CUBIC_SPLINE_INTERPOLATION;
}

//------------------------------------------------------------------------------

void ANEOSBirchMurnaghanDebyeModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 11, father);

  new ClassDouble<ANEOSBirchMurnaghanDebyeModelData>(ca, "DensityAtZeroKelvin", this, 
                       &ANEOSBirchMurnaghanDebyeModelData::zeroKelvinDensity);
  new ClassDouble<ANEOSBirchMurnaghanDebyeModelData>(ca, "ReferenceBulkModulus", this, 
                       &ANEOSBirchMurnaghanDebyeModelData::b0);
  new ClassDouble<ANEOSBirchMurnaghanDebyeModelData>(ca, "ReferenceBulkModulusDerivative", this, 
                       &ANEOSBirchMurnaghanDebyeModelData::b0prime);
  new ClassDouble<ANEOSBirchMurnaghanDebyeModelData>(ca, "EnergyShift", this, 
                       &ANEOSBirchMurnaghanDebyeModelData::delta_e);
  new ClassDouble<ANEOSBirchMurnaghanDebyeModelData>(ca, "MolarMass", this, 
                       &ANEOSBirchMurnaghanDebyeModelData::molar_mass);
  new ClassDouble<ANEOSBirchMurnaghanDebyeModelData>(ca, "ReferenceTemperature", this,
                       &ANEOSBirchMurnaghanDebyeModelData::T0);
  new ClassDouble<ANEOSBirchMurnaghanDebyeModelData>(ca, "ReferenceSpecificInternalEnergy", this, 
                       &ANEOSBirchMurnaghanDebyeModelData::e0);
  new ClassDouble<ANEOSBirchMurnaghanDebyeModelData>(ca, "ReferenceGamma", this, 
                       &ANEOSBirchMurnaghanDebyeModelData::Gamma0);
  new ClassDouble<ANEOSBirchMurnaghanDebyeModelData>(ca, "ReferenceDensity", this, 
                       &ANEOSBirchMurnaghanDebyeModelData::rho0);

  new ClassDouble<ANEOSBirchMurnaghanDebyeModelData>(ca, "BoltzmannConstant", this,
                       &ANEOSBirchMurnaghanDebyeModelData::boltzmann_constant);

  new ClassToken<ANEOSBirchMurnaghanDebyeModelData> (ca, "DebyeFunctionEvaluation", this,
        reinterpret_cast<int ANEOSBirchMurnaghanDebyeModelData::*>
            (&ANEOSBirchMurnaghanDebyeModelData::debye_evaluation), 
        2, "OnTheFly", 0, "CubicSplineInterpolation", 1);
}

//------------------------------------------------------------------------------

HyperelasticityModelData::HyperelasticityModelData()
{
  type = NONE;
  stress_option = DEVIATOR_ONLY;
  youngs_modulus = 0.0;
  poissons_ratio = 0.0;
  C01 = 0.0;
}

//------------------------------------------------------------------------------

void HyperelasticityModelData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassToken<HyperelasticityModelData>(ca, "Type", this,
           reinterpret_cast<int HyperelasticityModelData::*>(&HyperelasticityModelData::type), 5,
           "None",     HyperelasticityModelData::NONE,
           "SaintVenantKirchhoff", HyperelasticityModelData::SAINTVENANT_KIRCHHOFF,
           "ModifiedSaintVenantKirchhoff", HyperelasticityModelData::MODIFIED_SAINTVENANT_KIRCHHOFF,
           "NeoHookean", HyperelasticityModelData::NEO_HOOKEAN,
           "MooneyRivlin", HyperelasticityModelData::MOONEY_RIVLIN);

  new ClassToken<HyperelasticityModelData>(ca, "StressTensor", this,
           reinterpret_cast<int HyperelasticityModelData::*>(&HyperelasticityModelData::stress_option), 2,
           "Full",     HyperelasticityModelData::FULL,
           "DeviatorOnly", HyperelasticityModelData::DEVIATOR_ONLY);

  new ClassDouble<HyperelasticityModelData>(ca, "YoungsModulus", this,
                                            &HyperelasticityModelData::youngs_modulus);
  new ClassDouble<HyperelasticityModelData>(ca, "PoissonsRatio", this,
                                            &HyperelasticityModelData::poissons_ratio);
  new ClassDouble<HyperelasticityModelData>(ca, "C01", this, &HyperelasticityModelData::C01);
}

//------------------------------------------------------------------------------

MaterialModelData::MaterialModelData()
{

  eos = STIFFENED_GAS;
  rhomin = 0.0; // By default, density cannot be zero or negative
  pmin = -DBL_MAX;   // By default, no clipping
  rhomax = DBL_MAX;
  pmax = DBL_MAX;

  failsafe_density = 1.0e-10; //a small positive number

}

//------------------------------------------------------------------------------

Assigner *MaterialModelData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 16, nullAssigner);

  new ClassToken<MaterialModelData>(ca, "EquationOfState", this,
                                 reinterpret_cast<int MaterialModelData::*>(&MaterialModelData::eos), 7,
                                 "StiffenedGas", MaterialModelData::STIFFENED_GAS, 
                                 "NobleAbelStiffenedGas", MaterialModelData::NOBLE_ABEL_STIFFENED_GAS, 
                                 "MieGruneisen", MaterialModelData::MIE_GRUNEISEN,
                                 "ExtendedMieGruneisen", MaterialModelData::EXTENDED_MIE_GRUNEISEN,
                                 "Tillotson", MaterialModelData::TILLOTSON,
                                 "JonesWilkinsLee", MaterialModelData::JWL,
                                 "ANEOSBirchMurnaghanDebye", MaterialModelData::ANEOS_BIRCH_MURNAGHAN_DEBYE);
  new ClassDouble<MaterialModelData>(ca, "DensityCutOff", this, &MaterialModelData::rhomin);
  new ClassDouble<MaterialModelData>(ca, "PressureCutOff", this, &MaterialModelData::pmin);
  new ClassDouble<MaterialModelData>(ca, "DensityUpperLimit", this, &MaterialModelData::rhomax);
  new ClassDouble<MaterialModelData>(ca, "PressureUpperLimit", this, &MaterialModelData::pmax);

  new ClassDouble<MaterialModelData>(ca, "DensityPrescribedAtFailure", this, &MaterialModelData::failsafe_density);

  sgModel.setup("StiffenedGasModel", ca);
  nasgModel.setup("NobleAbelStiffenedGasModel", ca);
  mgModel.setup("MieGruneisenModel", ca);
  mgextModel.setup("ExtendedMieGruneisenModel", ca);
  tillotModel.setup("TillotsonModel", ca);
  jwlModel.setup("JonesWilkinsLeeModel", ca);
  abmdModel.setup("ANEOSBirchMurnaghanDebyeModel", ca);

  viscosity.setup("ViscosityModel", ca);
  
  heat_diffusion.setup("HeatDiffusionModel", ca);
  
  hyperelasticity.setup("HyperelasticityModel", ca);
  
  return ca;

};

//------------------------------------------------------------------------------

ViscosityModelData::ViscosityModelData() {

  type = NONE;

  dynamicViscosity   = 1.716e-5; //air, Pa.s
  bulkViscosity      = 0.0;

  sutherlandConstant = 111.0; //air, Kelvin
  sutherlandT0 = 273.0;  //air, Kelvin
  sutherlandMu0 = 1.716e-5;   //air, Pa.s

  Cav = 0.8;
  Cth = 0.05;
}

//------------------------------------------------------------------------------

void ViscosityModelData::setup(const char *name, ClassAssigner *father) {

  ClassAssigner *ca = new ClassAssigner(name, 8, father);

  new ClassToken<ViscosityModelData>(ca, "Type", this,
                                     reinterpret_cast<int ViscosityModelData::*>(&ViscosityModelData::type), 4,
                                     "None",        ViscosityModelData::NONE,
                                     "Constant",        ViscosityModelData::CONSTANT,
                                     "Sutherland",      ViscosityModelData::SUTHERLAND,
                                     "Artificial", ViscosityModelData::ARTIFICIAL_RODIONOV);

  new ClassDouble<ViscosityModelData>(ca, "ReferenceTemperature", this,
                                      &ViscosityModelData::sutherlandT0);
  new ClassDouble<ViscosityModelData>(ca, "ReferenceDynamicViscosity", this,
                                      &ViscosityModelData::sutherlandMu0);
  new ClassDouble<ViscosityModelData>(ca, "SutherlandConstant", this,
                                      &ViscosityModelData::sutherlandConstant);

  new ClassDouble<ViscosityModelData>(ca, "DynamicViscosity", this,
                                      &ViscosityModelData::dynamicViscosity);
  new ClassDouble<ViscosityModelData>(ca, "BulkViscosity", this,
                                      &ViscosityModelData::bulkViscosity);

  new ClassDouble<ViscosityModelData>(ca, "Cav", this,
                                      &ViscosityModelData::Cav);
  new ClassDouble<ViscosityModelData>(ca, "Cth", this,
                                      &ViscosityModelData::Cth);

}

//------------------------------------------------------------------------------

HeatDiffusionModelData::HeatDiffusionModelData()
{
  type = NONE;
  diffusivity = 0.0;
}

//------------------------------------------------------------------------------

void HeatDiffusionModelData::setup(const char *name, ClassAssigner *father) {

  ClassAssigner *ca = new ClassAssigner(name, 2, father);

  new ClassToken<HeatDiffusionModelData>(ca, "Type", this,
           reinterpret_cast<int HeatDiffusionModelData::*>(&HeatDiffusionModelData::type), 2,
           "None",     HeatDiffusionModelData::NONE,
           "Constant", HeatDiffusionModelData::CONSTANT);

  new ClassDouble<HeatDiffusionModelData>(ca, "Diffusivity", this, &HeatDiffusionModelData::diffusivity);
} 

//------------------------------------------------------------------------------

MaterialTransitionData::MaterialTransitionData()
{
  from_id = -1;
  to_id = -1;
  temperature_lowerbound = 0.0;
  temperature_upperbound = DBL_MAX;
  pressure_lowerbound = -DBL_MAX;
  pressure_upperbound = DBL_MAX;

  latent_heat = 0.0;
}

//------------------------------------------------------------------------------

Assigner *MaterialTransitionData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 7, nullAssigner);

  new ClassInt<MaterialTransitionData>(ca, "FromMaterialID", this, 
          &MaterialTransitionData::from_id);

  new ClassInt<MaterialTransitionData>(ca, "ToMaterialID", this, 
          &MaterialTransitionData::to_id);

  new ClassDouble<MaterialTransitionData>(ca, "TemperatureLowerbound", this, 
          &MaterialTransitionData::temperature_lowerbound);

  new ClassDouble<MaterialTransitionData>(ca, "TemperatureUpperbound", this, 
          &MaterialTransitionData::temperature_upperbound);

  new ClassDouble<MaterialTransitionData>(ca, "PressureLowerbound", this, 
          &MaterialTransitionData::pressure_lowerbound);

  new ClassDouble<MaterialTransitionData>(ca, "PressureUpperbound", this, 
          &MaterialTransitionData::pressure_upperbound);

  new ClassDouble<MaterialTransitionData>(ca, "LatentHeat", this, 
          &MaterialTransitionData::latent_heat);

  return ca;

}

//------------------------------------------------------------------------------

EquationsData::EquationsData()
{

}

//------------------------------------------------------------------------------

void EquationsData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 3, father); 

  materials.setup("Material", ca);

  transitions.setup("MaterialTransition", ca);

  dummy_state.setup("DummyState", ca);
}

//------------------------------------------------------------------------------

ReconstructionData::ReconstructionData() 
{
  type = LINEAR;
  limiter = NONE;
  slopeNearInterface = NONZERO;

  generalized_minmod_coeff = 1.2; //The generalized MC Limiter

  varType = PRIMITIVE;
}

//------------------------------------------------------------------------------

void ReconstructionData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 6, father); 

  new ClassToken<ReconstructionData>
    (ca, "Type", this,
     reinterpret_cast<int ReconstructionData::*>(&ReconstructionData::type), 2,
     "Constant", 0, "Linear", 1);

  new ClassToken<ReconstructionData>
    (ca, "SlopeNearInterface", this,
     reinterpret_cast<int ReconstructionData::*>(&ReconstructionData::slopeNearInterface), 2,
     "Zero", 0, "NonZero", 1);

  new ClassToken<ReconstructionData>
    (ca, "Limiter", this,
     reinterpret_cast<int ReconstructionData::*>(&ReconstructionData::limiter), 3,
     "None", 0, "GeneralizedMinMod", 1, "VanAlbada", 2);

  new ClassDouble<ReconstructionData>(ca, "GeneralizedMinModCoefficient", this, 
    &ReconstructionData::generalized_minmod_coeff);

  new ClassToken<ReconstructionData>
    (ca, "VariableType", this,
     reinterpret_cast<int ReconstructionData::*>(&ReconstructionData::varType), 4,
     "Primitive", 0, "Conservative", 1, "PrimitiveCharacteristic", 2,
     "ConservativeCharacteristic", 3);

  fixes.setup("Fixes", ca);
}

//------------------------------------------------------------------------------

SmoothingData::SmoothingData()
{
  type = NONE;
  iteration = 1;
  sigma_factor = 1.0;

  frequency = -100;
  frequency_dt = -1.0;

  conservation = OFF;
  conservation_tol = 0.1;
}

//------------------------------------------------------------------------------

void SmoothingData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 7, father); 

  new ClassToken<SmoothingData>
    (ca, "Type", this,
     reinterpret_cast<int SmoothingData::*>(&SmoothingData::type), 3,
     "None", 0, "Box", 1, "Gaussian", 2);

  new ClassInt<SmoothingData>(ca, "NumberOfIterations", this, &SmoothingData::iteration);

  new ClassDouble<SmoothingData>(ca, "SigmaFactor", this, &SmoothingData::sigma_factor);

  new ClassInt<SmoothingData>(ca, "Frequency", this, &SmoothingData::frequency);

  new ClassDouble<SmoothingData>(ca, "TimeInterval", this, &SmoothingData::frequency_dt);

  new ClassToken<SmoothingData>
    (ca, "EnforceConservation", this,
     reinterpret_cast<int SmoothingData::*>(&SmoothingData::conservation), 2,
     "Off", 0, "On", 1);

  new ClassDouble<SmoothingData>(ca, "ConservationTolerance", this, &SmoothingData::conservation_tol);
}

//------------------------------------------------------------------------------

FixData::FixData()
{ }

//------------------------------------------------------------------------------

void FixData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 5, father);
  sphereMap.setup("Sphere", ca);
  parallelepipedMap.setup("Parallelepiped", ca);
  spheroidMap.setup("Spheroid", ca);
  cylinderconeMap.setup("CylinderAndCone", ca);
  cylindersphereMap.setup("CylinderWithSphericalCaps", ca);
}

//------------------------------------------------------------------------------

SchemeData::SchemeData() 
{
  flux = HLLC;

  delta = 0.2; //the coefficient in Harten's entropy fix (for Roe flux)
}

//------------------------------------------------------------------------------

void SchemeData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner* ca;
  ca = new ClassAssigner(name, 4, father);

  new ClassToken<SchemeData>
    (ca, "Flux", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::flux), 4,
     "Roe", 0, "LocalLaxFriedrichs", 1, "HLLC", 2, "Godunov", 3);

  new ClassDouble<SchemeData>(ca, "EntropyFixCoefficient", this, &SchemeData::delta);

  rec.setup("Reconstruction", ca);

  smooth.setup("Smoothing", ca);

}

//------------------------------------------------------------------------------

BoundarySchemeData::BoundarySchemeData()
{
  type = STEGER_WARMING;
}

//------------------------------------------------------------------------------

void BoundarySchemeData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<BoundarySchemeData>(ca, "Type", this,
         reinterpret_cast<int BoundarySchemeData::*>(&BoundarySchemeData::type), 5,
         "StegerWarming", 0, "ConstantExtrapolation", 1, "LinearExtrapolation", 2,
         "Ghidaglia", 3, "ModifiedGhidaglia", 4);

}

//------------------------------------------------------------------------------

LevelSetReinitializationData::LevelSetReinitializationData()
{
  frequency = -1;
  frequency_dt = -1.0;
  maxIts = 30;
  cfl = 0.8;
  convergence_tolerance = 2.0e-4;
  firstLayerTreatment = FIXED;
}

//------------------------------------------------------------------------------

void LevelSetReinitializationData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassInt<LevelSetReinitializationData>(ca, "Frequency", this, 
          &LevelSetReinitializationData::frequency);

  new ClassDouble<LevelSetReinitializationData>(ca, "TimeInterval", this, 
          &LevelSetReinitializationData::frequency_dt);

  new ClassInt<LevelSetReinitializationData>(ca, "MaxIts", this, 
          &LevelSetReinitializationData::maxIts);

  new ClassDouble<LevelSetReinitializationData>(ca, "CFL", this, 
          &LevelSetReinitializationData::cfl);

  new ClassDouble<LevelSetReinitializationData>(ca, "ConvergenceTolerance", this, 
          &LevelSetReinitializationData::convergence_tolerance);

  new ClassToken<LevelSetReinitializationData>(ca, "FirstLayerTreatment", this,
     reinterpret_cast<int LevelSetReinitializationData::*>(&LevelSetReinitializationData::firstLayerTreatment), 5,
     "Fixed", 0, "ConstrainedMethod1", 1, "ConstrainedMethod2", 2,
     "IterativelyConstrainedMethod1", 3, "IterativelyConstrainedMethod2", 4);

}

//------------------------------------------------------------------------------

PrescribedMotionData::PrescribedMotionData()
{
  materialid = -1;

  velocity_x = 0.0;
  velocity_y = 0.0;
  velocity_z = 0.0;

  velocity_time_history = "";
}

//------------------------------------------------------------------------------

Assigner *PrescribedMotionData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 5, nullAssigner);

  new ClassInt<PrescribedMotionData>(ca, "MaterialID", this, 
          &PrescribedMotionData::materialid);

  new ClassDouble<PrescribedMotionData>(ca, "VelocityX", this, 
          &PrescribedMotionData::velocity_x);
  new ClassDouble<PrescribedMotionData>(ca, "VelocityY", this, 
          &PrescribedMotionData::velocity_y);
  new ClassDouble<PrescribedMotionData>(ca, "VelocityZ", this, 
          &PrescribedMotionData::velocity_z);

  new ClassStr<PrescribedMotionData>(ca, "VelocityTimeHistoryFile", this, 
          &PrescribedMotionData::velocity_time_history);

  return ca;
}

//------------------------------------------------------------------------------

LevelSetSchemeData::LevelSetSchemeData() 
{
  materialid = -1;

  flux = ROE;

  bandwidth = INT_MAX;

  solver = FINITE_DIFFERENCE;

  fd = UPWIND_CENTRAL_3;

  bc_x0   = ZERO_NEUMANN;
  bc_xmax = ZERO_NEUMANN;
  bc_y0   = ZERO_NEUMANN;
  bc_ymax = ZERO_NEUMANN;
  bc_z0   = ZERO_NEUMANN;
  bc_zmax = ZERO_NEUMANN;

  delta = 0.2; //the coefficient in Harten's entropy fix.

  init = DISTANCE_CALCULATION;
}

//------------------------------------------------------------------------------

Assigner *LevelSetSchemeData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 14, nullAssigner);

  new ClassInt<LevelSetSchemeData>(ca, "MaterialID", this, 
    &LevelSetSchemeData::materialid);

  new ClassToken<LevelSetSchemeData>
    (ca, "Solver", this,
     reinterpret_cast<int LevelSetSchemeData::*>(&LevelSetSchemeData::solver), 2,
     "FiniteVolume", 0, "FiniteDifference", 1);

  new ClassToken<LevelSetSchemeData>
    (ca, "Flux", this,
     reinterpret_cast<int LevelSetSchemeData::*>(&LevelSetSchemeData::flux), 3,
     "Roe", 0, "LocalLaxFriedrichs", 1, "Upwind", 2);

  new ClassDouble<LevelSetSchemeData>(ca, "EntropyFixCoefficient", this, &LevelSetSchemeData::delta);

  new ClassInt<LevelSetSchemeData>(ca, "Bandwidth", this, &LevelSetSchemeData::bandwidth);

  new ClassToken<LevelSetSchemeData>(ca, "BoundaryConditionX0", this,
          reinterpret_cast<int LevelSetSchemeData::*>(&LevelSetSchemeData::bc_x0), 4,
          "None", 0, "ZeroNeumann", 1, "LinearExtrapolation", 2, "NonNegative", 3);
  new ClassToken<LevelSetSchemeData>(ca, "BoundaryConditionXmax", this,
          reinterpret_cast<int LevelSetSchemeData::*>(&LevelSetSchemeData::bc_xmax), 4,
          "None", 0, "ZeroNeumann", 1, "LinearExtrapolation", 2, "NonNegative", 3);
  new ClassToken<LevelSetSchemeData>(ca, "BoundaryConditionY0", this,
          reinterpret_cast<int LevelSetSchemeData::*>(&LevelSetSchemeData::bc_y0), 4,
          "None", 0, "ZeroNeumann", 1, "LinearExtrapolation", 2, "NonNegative", 3);
  new ClassToken<LevelSetSchemeData>(ca, "BoundaryConditionYmax", this,
          reinterpret_cast<int LevelSetSchemeData::*>(&LevelSetSchemeData::bc_ymax), 4,
          "None", 0, "ZeroNeumann", 1, "LinearExtrapolation", 2, "NonNegative", 3);
  new ClassToken<LevelSetSchemeData>(ca, "BoundaryConditionZ0", this,
          reinterpret_cast<int LevelSetSchemeData::*>(&LevelSetSchemeData::bc_z0), 4,
          "None", 0, "ZeroNeumann", 1, "LinearExtrapolation", 2, "NonNegative", 3);
  new ClassToken<LevelSetSchemeData>(ca, "BoundaryConditionZmax", this,
          reinterpret_cast<int LevelSetSchemeData::*>(&LevelSetSchemeData::bc_zmax), 4,
          "None", 0, "ZeroNeumann", 1, "LinearExtrapolation", 2, "NonNegative", 3);

  new ClassToken<LevelSetSchemeData>(ca, "Initialization", this,
          reinterpret_cast<int LevelSetSchemeData::*>(&LevelSetSchemeData::init), 2,
          "DistanceCalculation", 0, "Reinitialization", 1);

  rec.setup("Reconstruction", ca);

  reinit.setup("Reinitialization", ca);

  return ca;
}

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------

SchemesData::SchemesData() 
{

}

//------------------------------------------------------------------------------

void SchemesData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  ns.setup("NavierStokes", ca);

  bc.setup("Boundaries", ca);

  ls.setup("LevelSet", ca);

  pm.setup("PrescribedMotion", ca);

}

//------------------------------------------------------------------------------

ExactRiemannSolverData::ExactRiemannSolverData()
{
  maxIts_main = 200;
  maxIts_bracket = 100;
  maxIts_shock = 200;
  numSteps_rarefaction = 200;
  tol_main = 1.0e-4; //applied to both pressure and velocity
  tol_shock = 1.0e-12; //a density tolerance (non-D)
  tol_rarefaction = 1.0e-8; //a density and velocity tolerance (non-D); 
  min_pressure = -1.0e8;
  failure_threshold = 0.2;
  pressure_at_failure = 1.0e-8;

  // Experimental
  surface_tension = NO;
  surface_tension_coefficient = 0.;
  surface_tension_materialid = 1;
}

//------------------------------------------------------------------------------

void ExactRiemannSolverData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 13, father);

  new ClassInt<ExactRiemannSolverData>(ca, "MaxIts", this, 
                                       &ExactRiemannSolverData::maxIts_main);

  new ClassInt<ExactRiemannSolverData>(ca, "MaxItsBracketing", this, 
                                       &ExactRiemannSolverData::maxIts_bracket);

  new ClassInt<ExactRiemannSolverData>(ca, "MaxItsShock", this, 
                                       &ExactRiemannSolverData::maxIts_shock);

  new ClassInt<ExactRiemannSolverData>(ca, "IntegrationSteps", this, 
                                       &ExactRiemannSolverData::numSteps_rarefaction);

  new ClassDouble<ExactRiemannSolverData>(ca, "Tolerance", this, 
                                          &ExactRiemannSolverData::tol_main);

  new ClassDouble<ExactRiemannSolverData>(ca, "ToleranceShock", this, 
                                          &ExactRiemannSolverData::tol_shock);

  new ClassDouble<ExactRiemannSolverData>(ca, "ToleranceRarefaction", this, 
                                          &ExactRiemannSolverData::tol_rarefaction);

  new ClassDouble<ExactRiemannSolverData>(ca, "MinPressure", this,
                                          &ExactRiemannSolverData::min_pressure);

  new ClassDouble<ExactRiemannSolverData>(ca, "FailureThreshold", this,
                                          &ExactRiemannSolverData::failure_threshold);

  new ClassDouble<ExactRiemannSolverData>(ca, "PrescribedPressureUponFailure", this,
                                          &ExactRiemannSolverData::pressure_at_failure);

  // Experimental 
  
  new ClassToken<ExactRiemannSolverData>(ca, "SurfaceTension", this,
                                         reinterpret_cast<int ExactRiemannSolverData::*>
                                         (&ExactRiemannSolverData::surface_tension), 2,
                                         "No", 0, "Yes", 1);

  new ClassDouble<ExactRiemannSolverData>(ca, "SurfaceTensionCoefficient", this,
                                          &ExactRiemannSolverData::surface_tension_coefficient);

  new ClassInt<ExactRiemannSolverData>(ca, "SurfaceTensionMaterialID", this,
                                       &ExactRiemannSolverData::surface_tension_materialid);
}

//------------------------------------------------------------------------------

MultiPhaseData::MultiPhaseData()
{
  flux = NUMERICAL;
  recon = CONSTANT;

  phasechange_type = RIEMANN_SOLUTION;
  phasechange_dir  = UPWIND;

  riemann_normal = MESH;

  latent_heat_transfer = Off;

  levelset_correction_frequency = -1;

  apply_failsafe_density = On;

  conRec_depth = 0.0;
}

//------------------------------------------------------------------------------

void MultiPhaseData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 9, father);

  new ClassToken<MultiPhaseData>
    (ca, "Flux", this,
     reinterpret_cast<int MultiPhaseData::*>(&MultiPhaseData::flux), 3,
     "Exact", 0, "Numerical", 1, "LocalLaxFriedrichs", 2);

  new ClassToken<MultiPhaseData>
    (ca, "ReconstructionAtInterface", this,
     reinterpret_cast<int MultiPhaseData::*>(&MultiPhaseData::recon), 2,
     "Constant", 0, "Linear", 1);

  new ClassToken<MultiPhaseData>
    (ca, "PhaseChange", this,
     reinterpret_cast<int MultiPhaseData::*>(&MultiPhaseData::phasechange_type), 2,
     "RiemannSolution", 0, "Extrapolation", 1);

  new ClassToken<MultiPhaseData>
    (ca, "PhaseChangeDirection", this,
     reinterpret_cast<int MultiPhaseData::*>(&MultiPhaseData::phasechange_dir), 2,
     "All", 0, "Upwind", 1);

  new ClassToken<MultiPhaseData>
    (ca, "RiemannNormal", this,
     reinterpret_cast<int MultiPhaseData::*>(&MultiPhaseData::riemann_normal), 3,
     "LevelSet", 0, "Mesh", 1, "Average", 2);

  new ClassToken<MultiPhaseData>
    (ca, "LatentHeatTransfer", this,
     reinterpret_cast<int MultiPhaseData::*>(&MultiPhaseData::latent_heat_transfer), 2,
     "Off", 0, "On", 1);

  new ClassInt<MultiPhaseData>(ca, "LevelSetCorrectionFrequency", 
        this, &MultiPhaseData::levelset_correction_frequency);


  new ClassToken<MultiPhaseData>
    (ca, "ApplyFailSafeDensity", this,
     reinterpret_cast<int MultiPhaseData::*>(&MultiPhaseData::apply_failsafe_density), 2,
     "Off", 0, "On", 1);

  new ClassDouble<MultiPhaseData>(ca, "ConstantReconstructionDepth", 
        this, &MultiPhaseData::conRec_depth);

}

//------------------------------------------------------------------------------

ExplicitData::ExplicitData()
{
  type = RUNGE_KUTTA_2;
}

//------------------------------------------------------------------------------

void ExplicitData::setup(const char *name, ClassAssigner *father)
{

 ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<ExplicitData>
    (ca, "Type", this,
     reinterpret_cast<int ExplicitData::*>(&ExplicitData::type), 3,
     "ForwardEuler", 0, "RungeKutta2", 1, "RungeKutta3", 2);

}

//------------------------------------------------------------------------------

TsData::TsData()
{

  type = EXPLICIT;
  maxIts = INT_MAX;
  timestep = -1.0;
  cfl = 0.5;
  maxTime = 1e10;

  convergence_tolerance = -1.0; //!< activated only for steady-state computations
  local_dt = NO;

}

//------------------------------------------------------------------------------

void TsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 8, father);

  new ClassToken<TsData>(ca, "Type", this,
                         reinterpret_cast<int TsData::*>(&TsData::type), 2,
                         "Explicit", 0, "Implicit", 1);
  new ClassInt<TsData>(ca, "MaxIts", this, &TsData::maxIts);
  new ClassDouble<TsData>(ca, "TimeStep", this, &TsData::timestep);
  new ClassDouble<TsData>(ca, "CFL", this, &TsData::cfl);
  new ClassDouble<TsData>(ca, "MaxTime", this, &TsData::maxTime);

  new ClassDouble<TsData>(ca, "ConvergenceTolerance", this, &TsData::convergence_tolerance);
  new ClassToken<TsData>(ca, "LocalTimeStepping", this,
                         reinterpret_cast<int TsData::*>(&TsData::local_dt), 2,
                         "Off", 0, "On", 1);


  expl.setup("Explicit", ca);
}

//------------------------------------------------------------------------------

RectangleData::RectangleData()
{
  
  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  normal_x = 0.0;
  normal_y = 0.0;
  normal_z = 0.0;
  a = -1.0;
  b = -1.0;

}

//------------------------------------------------------------------------------

Assigner *RectangleData::getAssigner()
{
  
  ClassAssigner *ca = new ClassAssigner("normal", 9, nullAssigner);
  
  new ClassDouble<RectangleData> (ca, "Center_x", this, &RectangleData::cen_x);
  new ClassDouble<RectangleData> (ca, "Center_y", this, &RectangleData::cen_y);
  new ClassDouble<RectangleData> (ca, "Center_z", this, &RectangleData::cen_z);
  new ClassDouble<RectangleData> (ca, "Normal_x", this, &RectangleData::normal_x);
  new ClassDouble<RectangleData> (ca, "Normal_y", this, &RectangleData::normal_y);
  new ClassDouble<RectangleData> (ca, "Normal_z", this, &RectangleData::normal_z);
  new ClassDouble<RectangleData> (ca, "Dimension1", this, &RectangleData::a);
  new ClassDouble<RectangleData> (ca, "Dimension2", this, &RectangleData::b);
  
  state.setup("BoundaryState", ca);
  
  return ca;
}

//------------------------------------------------------------------------------


DiskData::DiskData()
{
  
  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  normal_x = 0.0;
  normal_y = 0.0;
  normal_z = 0.0;
  radius = -1.0;

}

//------------------------------------------------------------------------------

Assigner *DiskData::getAssigner()
{
  
  ClassAssigner *ca = new ClassAssigner("normal", 8, nullAssigner);
  
  new ClassDouble<DiskData> (ca, "Center_x", this, &DiskData::cen_x);
  new ClassDouble<DiskData> (ca, "Center_y", this, &DiskData::cen_y);
  new ClassDouble<DiskData> (ca, "Center_z", this, &DiskData::cen_z);
  new ClassDouble<DiskData> (ca, "Normal_x", this, &DiskData::normal_x);
  new ClassDouble<DiskData> (ca, "Normal_y", this, &DiskData::normal_y);
  new ClassDouble<DiskData> (ca, "Normal_z", this, &DiskData::normal_z);
  new ClassDouble<DiskData> (ca, "Radius", this, &DiskData::radius);
  
  state.setup("BoundaryState", ca);
  
  return ca;
}

//------------------------------------------------------------------------------

void MultiBoundaryConditionsData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 2, father);
  diskMap.setup("Disk", ca);
  rectangleMap.setup("Rectangle", ca);
}

//------------------------------------------------------------------------------

BcsData::BcsData()
{

}

//------------------------------------------------------------------------------

void BcsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 7, father);

  inlet.setup("Inlet", ca);
  //Inside the code, Farfield0 = Farfield = Inlet, Farfield1 = Outlet
  inlet.setup("Farfield0", ca);
  inlet.setup("Farfield", ca);
  outlet.setup("Outlet", ca);
  outlet.setup("Farfield1", ca);
  wall.setup("Wall", ca);

  multiBoundaryConditions.setup("GeometricEntities2D");

}

//------------------------------------------------------------------------------

BcsWallData::BcsWallData()
{

  type = ADIABATIC;
  temperature = -1.0;

}

//------------------------------------------------------------------------------

void BcsWallData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 2, father);

  new ClassToken<BcsWallData>(ca, "Type", this,
                              reinterpret_cast<int BcsWallData::*>(&BcsWallData::type), 2,
                              "Isothermal", 0, "Adiabatic", 1);
  new ClassDouble<BcsWallData>(ca, "Temperature", this, &BcsWallData::temperature);

}

//------------------------------------------------------------------------------

FloodIcData::FloodIcData()
{
  source_x = DBL_MAX; //DBL_MAX means the user did not specify it.
  source_y = DBL_MAX;
  source_z = DBL_MAX;

  waterline_x = 0.0;
  waterline_y = 0.0;
  waterline_z = 0.0;

  gx = 0.0;
  gy = 0.0;
  gz = 0.0;
}

//------------------------------------------------------------------------------

void FloodIcData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 10, father);

  new ClassDouble<FloodIcData>(ca, "Source_x", this, &FloodIcData::source_x);
  new ClassDouble<FloodIcData>(ca, "Source_y", this, &FloodIcData::source_y);
  new ClassDouble<FloodIcData>(ca, "Source_z", this, &FloodIcData::source_z);

  new ClassDouble<FloodIcData>(ca, "Waterline_x", this, &FloodIcData::waterline_x);
  new ClassDouble<FloodIcData>(ca, "Waterline_y", this, &FloodIcData::waterline_y);
  new ClassDouble<FloodIcData>(ca, "Waterline_z", this, &FloodIcData::waterline_z);

  new ClassDouble<FloodIcData>(ca, "Gravity_x", this, &FloodIcData::gx);
  new ClassDouble<FloodIcData>(ca, "Gravity_y", this, &FloodIcData::gy);
  new ClassDouble<FloodIcData>(ca, "Gravity_z", this, &FloodIcData::gz);

  waterline_ic.setup("InitialState");
}

//------------------------------------------------------------------------------

IcData::IcData()
{
  user_specified_ic = "";

  apply_user_file_before_geometries = NO;

  rbf = MULTIQUADRIC;

  type = NONE;

  for(int i=0; i<SIZE; i++)
    specified[i] = 0;
}

//------------------------------------------------------------------------------

void IcData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassStr<IcData>(ca, "UserDataFile", this, &IcData::user_specified_ic);

  new ClassToken<IcData> (ca, "ApplyUserDataBeforeGeometricEntities", this,
        reinterpret_cast<int IcData::*>(&IcData::apply_user_file_before_geometries), 2,
        "No", 0, "Yes", 1);

  new ClassToken<IcData> (ca, "InterpolationFunction", this,
        reinterpret_cast<int IcData::*>(&IcData::rbf), 4, "Multiquadric", 0, 
        "InverseMultiquadric", 1, "ThinPlateSpline", 2, "Gaussian", 3);

  default_ic.setup("DefaultInitialState");

  multiInitialConditions.setup("GeometricEntities");

  floodIc.setup("Flood");

}

//------------------------------------------------------------------------------

void IcData::readUserSpecifiedIC()
{
  if(!strcmp(user_specified_ic, "")) //no user_specified_ic
    return;

  std::fstream input;
  input.open(user_specified_ic, std::fstream::in);
  if (!input.is_open()) {
    print_error("*** Error: Cannot open user-specified initial condition file %s.\n", user_specified_ic);
    exit_mpi();
  } else
    print("- Reading user-specified initial condition file: %s.\n", user_specified_ic);

  std::string word, line;

  input.ignore(512,'\n'); //done with line 1 (not read, user's comment)

  // Read the second line of user-specified file
  input.ignore(2,' '); //This line should start with ## (so the file can be plotted
                       //directly using gnuplot). It must then specify the type of initial cond.

  input >> word;
  // just comparing the first 4 letters is sufficient
  if(!(word.compare(0,4,"Planar",0,4) && 
       word.compare(0,4,"PLANAR",0,4) && 
       word.compare(0,4,"planar",0,4))) {
    type = PLANAR;
    input.ignore(256,'\n'); //done with line 2
    readUserSpecifiedIC_Planar(input);
  } 
  else if(!(word.compare(0,4,"Cylindrical",0,4) && 
            word.compare(0,4,"CYLINDRICAL",0,4) && 
            word.compare(0,4,"cylindrical",0,4))) {
    type = CYLINDRICAL;
    input.ignore(256,'\n'); //done with line 2
    readUserSpecifiedIC_Cylindrical(input);
  } 
  else if(!(word.compare(0,4,"Spherical",0,4) && 
            word.compare(0,4,"SPHERICAL",0,4) && 
            word.compare(0,4,"spherical",0,4))) {
    type = SPHERICAL;
    input.ignore(256,'\n'); //done with line 2
    readUserSpecifiedIC_Spherical(input);
  } 
  else if(!(word.compare(0,4,"GeneralCylindrical",0,4) && 
            word.compare(0,4,"GENERALCYLINDRICAL",0,4) && 
            word.compare(0,4,"generalcylindrical",0,4))) {
    type = GENERALCYLINDRICAL;
    input.ignore(256,'\n'); //done with line 2
    readUserSpecifiedIC_GeneralCylindrical(input);
  } 
  else {
    print_error("*** Error: Unknown initial condition type %s.\n", word);
    exit_mpi();
  }

  input.close();

}

//------------------------------------------------------------------------------

// Read planar initial condition from file
void IcData::readUserSpecifiedIC_Planar(std::fstream &input)
{
  std::string word, line;

  // Read the third line of user-specified file
  input.ignore(2,' '); //This line should start with ## 
                       //It must then contain 3 real numbers corresponding to the (x,y,z) coordinates
                       //of the "0" in this data file within the actual mesh
  input >> x0[0] >> x0[1] >> x0[2];
  input.ignore(256,'\n'); //done with this line.

  // Read the next line
  input.ignore(2,' '); //This line should start with ##
                       //It must contain 3 real numbers corresponding to the direction of "+x" in the 1D data
  input >> dir[0] >> dir[1] >> dir[2];
  dir /= dir.norm();
  input.ignore(256,'\n'); //done with this line

  // Read the next line, which specifies the variable of each column
  input.ignore(2,' '); //This line should start with ##
                       //It must contain between 2 and SIZE words indicating the field
  for(int i=0; i<SIZE; i++)
    specified[i] = 0;
  std::map<int,int> column2var; //maps the column number to "Vars" index 
  getline(input, line);

  std::istringstream iss(line);
  int column = 0; 
  for(;iss>>word;) {
    // just check the first four letters
    if(!(word.compare(0,4,"Coordinate",0,4) && 
         word.compare(0,4,"COORDINATE",0,4) && 
         word.compare(0,4,"coordinate",0,4))) {
      column2var[column] = COORDINATE;
      specified[COORDINATE] = 1;
    } 
    else if(!(word.compare(0,4,"Density",0,4) && 
              word.compare(0,4,"DENSITY",0,4) && 
              word.compare(0,4,"density",0,4))) {
      column2var[column] = DENSITY;
      specified[DENSITY] = 1;
    } 
    else if(!(word.compare(0,4,"Velocity",0,4) && 
              word.compare(0,4,"VELOCITY",0,4) && 
              word.compare(0,4,"velocity",0,4))) {
      column2var[column] = VELOCITY;
      specified[VELOCITY] = 1;
    } 
    else if(!(word.compare(0,4,"Pressure",0,4) && 
              word.compare(0,4,"PRESSURE",0,4) && 
              word.compare(0,4,"pressure",0,4))) {
      column2var[column] = PRESSURE;
      specified[PRESSURE] = 1;
    } 
    else if(!(word.compare(0,4,"LevelSet",0,4) && 
              word.compare(0,4,"LEVELSET",0,4) && 
              word.compare(0,4,"levelset",0,4))) {
      column2var[column] = LEVELSET;
      specified[LEVELSET] = 1;
    } 
    else if(!(word.compare(0,4,"MaterialID",0,4) && 
              word.compare(0,4,"MATERIALID",0,4) && 
              word.compare(0,4,"materialid",0,4))) {
      column2var[column] = MATERIALID;
      specified[MATERIALID] = 1;
    } 
    else if(!(word.compare(0,4,"Temperature",0,4) && 
              word.compare(0,4,"TEMPERATURE",0,4) && 
              word.compare(0,4,"temperature",0,4))) {
      column2var[column] = TEMPERATURE;
      specified[TEMPERATURE] = 1;
    } else {
      print_error("*** Error: I do not understand the word '%s' in the user-specified initial condition file.\n", word.c_str());
      exit_mpi();
    }
    column++;
  }

  if(column<2 || !specified[COORDINATE]) {
    print_error("*** Error: Need additional data in the initial condition file.\n");
    exit_mpi();
  }

  // Now start reading the actual data (until end of file)
  double data;
  for(int r=0; r<INT_MAX; r++) {
    getline(input, line);
    std::istringstream is(line);
    for(int i=0; i<column; i++) {
      is >> data; 
      if(is.fail())
        break;
      //cout << "row " << r << ", column " << i << ": " << data << endl;
      user_data[column2var[i]].push_back(data);
    }
    if(input.eof())
      break;
  }

}

//------------------------------------------------------------------------------

// Read cylindrical initial condition from file
void IcData::readUserSpecifiedIC_Cylindrical(std::fstream &input)
{
  std::string word, line;

  // Read the third line of user-specified file
  input.ignore(2,' '); //This line should start with ## 
                       //It must then contain 3 real numbers corresponding to the (x,y,z) coordinates
                       //of the "0" in this data file within the actual mesh
  input >> x0[0] >> x0[1] >> x0[2];
  input.ignore(256,'\n'); //done with line 2

  // Read the next line
  input.ignore(2,' '); //This line should start with ##
                       //It must contain 3 real numbers corresponding to the positive direction of axis of symmetry
  input >> dir[0] >> dir[1] >> dir[2];
  dir /= dir.norm();
  input.ignore(256,'\n'); //done with this line

  // Read the next line, which specifies the variable of each column
  input.ignore(2,' '); //This line should start with ##
                       //It must contain between 2 and SIZE words indicating the field
  for(int i=0; i<SIZE; i++)
    specified[i] = 0;
  std::map<int,int> column2var; //maps the column number to "Vars" index 
  getline(input, line);

  std::istringstream iss(line);
  int column = 0; 
  for(;iss>>word;) {
    // just check the first four letters
    if(!(word.compare(0,4,"Coordinate",0,4) && 
         word.compare(0,4,"COORDINATE",0,4) && 
         word.compare(0,4,"coordinate",0,4))) {
      column2var[column] = COORDINATE;
      specified[COORDINATE] = 1;
    } 
    else if(!(word.compare(0,4,"Density",0,4) && 
              word.compare(0,4,"DENSITY",0,4) && 
              word.compare(0,4,"density",0,4))) {
      column2var[column] = DENSITY;
      specified[DENSITY] = 1;
    } 
    else if(!(word.compare(0,4,"Velocity",0,4) && 
              word.compare(0,4,"VELOCITY",0,4) && 
              word.compare(0,4,"velocity",0,4))) {
      column2var[column] = VELOCITY;
      specified[VELOCITY] = 1;
    } 
    else if(!(word.compare(0,4,"Pressure",0,4) && 
              word.compare(0,4,"PRESSURE",0,4) && 
              word.compare(0,4,"pressure",0,4))) {
      column2var[column] = PRESSURE;
      specified[PRESSURE] = 1;
    } 
    else if(!(word.compare(0,4,"LevelSet",0,4) && 
              word.compare(0,4,"LEVELSET",0,4) && 
              word.compare(0,4,"levelset",0,4))) {
      column2var[column] = LEVELSET;
      specified[LEVELSET] = 1;
    } 
    else if(!(word.compare(0,4,"MaterialID",0,4) && 
              word.compare(0,4,"MATERIALID",0,4) && 
              word.compare(0,4,"materialid",0,4))) {
      column2var[column] = MATERIALID;
      specified[MATERIALID] = 1;
    } 
    else if(!(word.compare(0,4,"Temperature",0,4) && 
              word.compare(0,4,"TEMPERATURE",0,4) && 
              word.compare(0,4,"temperature",0,4))) {
      column2var[column] = TEMPERATURE;
      specified[TEMPERATURE] = 1;
    } else {
      print_error("*** Error: I do not understand the word '%s' in the user-specified initial condition file.\n", word.c_str());
      exit_mpi();
    }
    column++;
  }

  if(column<2 || !specified[COORDINATE]) {
    print_error("*** Error: Need additional data in the initial condition file.\n");
    exit_mpi();
  }


  // Read the next line: should start with ##, then say "AXIAL"
  input.ignore(2,' '); 
  input >> word;
  if(word.compare(0,4,"Axial",0,4) && 
     word.compare(0,4,"AXIAL",0,4) && 
     word.compare(0,4,"axial",0,4)) {
    print_error("*** Error: Expect keyword 'Axial' in user-specified initial condition.\n");
    exit_mpi();
  } 
  input.ignore(256,'\n'); //done with this line


  // Now start reading the data in the axial direction
  double data;
  bool found_radial = false;
  for(int r=0; r<INT_MAX; r++) {

    getline(input, line);

    // check if we got the line that says "## Radial"
    std::istringstream check(line);
    check >> word; check >> word;
    if(!(word.compare(0,4,"Radial",0,4) && 
         word.compare(0,4,"RADIAL",0,4) && 
         word.compare(0,4,"radial",0,4))) {
      found_radial = true;
      break; 
    }

    std::istringstream is(line);
    for(int i=0; i<column; i++) {
      is >> data; 
      if(is.fail())
        break;
      //cout << "row " << r << ", column " << i << ": " << data << endl;
      user_data[column2var[i]].push_back(data);
    }

    if(input.eof())
      break;
  }

  if(found_radial) { //read radial variation
    for(int r=0; r<INT_MAX; r++) {
      getline(input, line);

      std::istringstream is(line);
      for(int i=0; i<column; i++) {
        is >> data; 
        if(is.fail())
          break;
        //cout << "row " << r << ", column " << i << ": " << data << endl;
        user_data2[column2var[i]].push_back(data);
      }

      if(input.eof())
        break;
    }
  }

}


//------------------------------------------------------------------------------

// Read general cylindrical initial condition (2D) from file
void IcData::readUserSpecifiedIC_GeneralCylindrical(std::fstream &input)
{
  std::string word, line;

  // Read the third line of user-specified file
  input.ignore(2,' '); //This line should start with ## 
                       //It must then contain 3 real numbers corresponding to the (x,y,z) coordinates
                       //of the "0" in this data file within the actual mesh
  input >> x0[0] >> x0[1] >> x0[2];
  input.ignore(256,'\n'); //done with line 2

  // Read the next line
  input.ignore(2,' '); //This line should start with ##
                       //It must contain 3 real numbers corresponding to the positive direction of axis of symmetry
  input >> dir[0] >> dir[1] >> dir[2];
  dir /= dir.norm();
  input.ignore(256,'\n'); //done with this line

  // Read the next line
  input.ignore(2,' '); //This line should start with ##
                       //It must contain 4 real numbers corresponding to xmin, xmax, rmin, rmax --- the bounding box
                       //for interpolation
  input >> xmin[0] >> xmax[0] >> xmin[1] >> xmax[1];
  input.ignore(256,'\n'); //done with this line

  // Read the next line, which specifies the variable of each column
  input.ignore(2,' '); //This line should start with ##
                       //It must contain between 2 and SIZE words indicating the field
  for(int i=0; i<SIZE; i++)
    specified[i] = 0;
  std::map<int,int> column2var; //maps the column number to "Vars" index 
  getline(input, line);

  std::istringstream iss(line);
  int column = 0; 
  for(;iss>>word;) {
    // just check the first four letters
    if(!(word.compare(0,10,"AxialCoordinate",0,10) && 
         word.compare(0,10,"AXIALCOORDINATE",0,10) && 
         word.compare(0,10,"axialcoordinate",0,10))) {
      column2var[column] = COORDINATE;
      specified[COORDINATE] = 1;
    } 
    else if(!(word.compare(0,10,"RadialCoordinate",0,10) && 
         word.compare(0,10,"RADIALCOORDINATE",0,10) && 
         word.compare(0,10,"radialcoordinate",0,10))) {
      column2var[column] = RADIALCOORDINATE;
      specified[RADIALCOORDINATE] = 1;
    } 
    else if(!(word.compare(0,4,"Density",0,4) && 
              word.compare(0,4,"DENSITY",0,4) && 
              word.compare(0,4,"density",0,4))) {
      column2var[column] = DENSITY;
      specified[DENSITY] = 1;
    } 
    else if(!(word.compare(0,10,"AxialVelocity",0,10) && 
              word.compare(0,10,"AXIALVELOCITY",0,10) && 
              word.compare(0,10,"axialvelocity",0,10))) {
      column2var[column] = VELOCITY;
      specified[VELOCITY] = 1;
    } 
    else if(!(word.compare(0,10,"RadialVelocity",0,10) && 
              word.compare(0,10,"RADIALVELOCITY",0,10) && 
              word.compare(0,10,"radialvelocity",0,10))) {
      column2var[column] = RADIALVELOCITY;
      specified[RADIALVELOCITY] = 1;
    } 
    else if(!(word.compare(0,4,"Pressure",0,4) && 
              word.compare(0,4,"PRESSURE",0,4) && 
              word.compare(0,4,"pressure",0,4))) {
      column2var[column] = PRESSURE;
      specified[PRESSURE] = 1;
    } 
    else if(!(word.compare(0,4,"LevelSet",0,4) && 
              word.compare(0,4,"LEVELSET",0,4) && 
              word.compare(0,4,"levelset",0,4))) {
      column2var[column] = LEVELSET;
      specified[LEVELSET] = 1;
    } 
    else if(!(word.compare(0,4,"MaterialID",0,4) && 
              word.compare(0,4,"MATERIALID",0,4) && 
              word.compare(0,4,"materialid",0,4))) {
      column2var[column] = MATERIALID;
      specified[MATERIALID] = 1;
    } 
    else if(!(word.compare(0,4,"Temperature",0,4) && 
              word.compare(0,4,"TEMPERATURE",0,4) && 
              word.compare(0,4,"temperature",0,4))) {
      column2var[column] = TEMPERATURE;
      specified[TEMPERATURE] = 1;
    } else {
      print_error("*** Error: I do not understand the word '%s' in the user-specified initial condition file.\n", word.c_str());
      exit_mpi();
    }
    column++;
  }

  if(column<3 || !specified[COORDINATE] || !specified[RADIALCOORDINATE]) {
    print_error("*** Error: Need additional data in the initial condition file.\n");
    exit_mpi();
  }

  // Now start reading the data in the axial direction
  double data;
  //bool found_radial = false;
  for(int r=0; r<INT_MAX; r++) {

    getline(input, line);

    std::istringstream is(line);
    for(int i=0; i<column; i++) {
      is >> data; 
      if(is.fail())
        break;
      //cout << "row " << r << ", column " << i << ": " << data << endl;
      user_data[column2var[i]].push_back(data);
    }

    if(input.eof())
      break;
  }

}


//------------------------------------------------------------------------------

// Read spherical initial condition from file
// Note: This is almost identical to reading planar data, except the variable "dir" does not need be read
void IcData::readUserSpecifiedIC_Spherical(std::fstream &input)
{
  std::string word, line;

  // Read the third line of user-specified file
  input.ignore(2,' '); //This line should start with ## 
                       //It must then contain 3 real numbers corresponding to the (x,y,z) coordinates
                       //of the "0" in this data file within the actual mesh
  input >> x0[0] >> x0[1] >> x0[2];
  input.ignore(256,'\n'); //done with line 2

  // Read the next line, which specifies the variable of each column
  input.ignore(2,' '); //This line should start with ##
                       //It must contain between 2 and SIZE words indicating the field
  for(int i=0; i<SIZE; i++)
    specified[i] = 0;
  std::map<int,int> column2var; //maps the column number to "Vars" index 
  getline(input, line);

  std::istringstream iss(line);
  int column = 0; 
  for(;iss>>word;) {
    // just check the first four letters
    if(!(word.compare(0,4,"Coordinate",0,4) && 
         word.compare(0,4,"COORDINATE",0,4) && 
         word.compare(0,4,"coordinate",0,4))) {
      column2var[column] = COORDINATE;
      specified[COORDINATE] = 1;
    } 
    else if(!(word.compare(0,4,"Density",0,4) && 
              word.compare(0,4,"DENSITY",0,4) && 
              word.compare(0,4,"density",0,4))) {
      column2var[column] = DENSITY;
      specified[DENSITY] = 1;
    } 
    else if(!(word.compare(0,4,"Velocity",0,4) && 
              word.compare(0,4,"VELOCITY",0,4) && 
              word.compare(0,4,"velocity",0,4))) {
      column2var[column] = VELOCITY;
      specified[VELOCITY] = 1;
    } 
    else if(!(word.compare(0,4,"Pressure",0,4) && 
              word.compare(0,4,"PRESSURE",0,4) && 
              word.compare(0,4,"pressure",0,4))) {
      column2var[column] = PRESSURE;
      specified[PRESSURE] = 1;
    } 
    else if(!(word.compare(0,4,"LevelSet",0,4) && 
              word.compare(0,4,"LEVELSET",0,4) && 
              word.compare(0,4,"levelset",0,4))) {
      column2var[column] = LEVELSET;
      specified[LEVELSET] = 1;
    } 
    else if(!(word.compare(0,4,"MaterialID",0,4) && 
              word.compare(0,4,"MATERIALID",0,4) && 
              word.compare(0,4,"materialid",0,4))) {
      column2var[column] = MATERIALID;
      specified[MATERIALID] = 1;
    } 
    else if(!(word.compare(0,4,"Temperature",0,4) && 
              word.compare(0,4,"TEMPERATURE",0,4) && 
              word.compare(0,4,"temperature",0,4))) {
      column2var[column] = TEMPERATURE;
      specified[TEMPERATURE] = 1;
    } else {
      print_error("*** Error: I do not understand the word '%s' in the user-specified initial condition file.\n", word.c_str());
      exit_mpi();
    }
    column++;
  }

  if(column<2 || !specified[COORDINATE]) {
    print_error("*** Error: Need additional data in the initial condition file.\n");
    exit_mpi();
  }

  // Now start reading the actual data (until end of file)
  int MaxRows = 200000; //should be a lot more than enough
  double data;
  for(int r=0; r<MaxRows; r++) {
    getline(input, line);
    std::istringstream is(line);
    for(int i=0; i<column; i++) {
      is >> data; 
      if(is.fail())
        break;
      //cout << "row " << r << ", column " << i << ": " << data << endl;
      user_data[column2var[i]].push_back(data);
    }
    if(input.eof())
      break;
  }

}

//------------------------------------------------------------------------------

LaserAbsorptionCoefficient::LaserAbsorptionCoefficient()
{
  materialid = 0;
  slope  = 0.0;
  T0     = 0.0;
  alpha0 = 0.0;
}

//------------------------------------------------------------------------------

Assigner* LaserAbsorptionCoefficient::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 4, nullAssigner);
  new ClassInt<LaserAbsorptionCoefficient>(ca, "MaterialID", this, &LaserAbsorptionCoefficient::materialid);
  new ClassDouble<LaserAbsorptionCoefficient>(ca, "Slope", this, &LaserAbsorptionCoefficient::slope);
  new ClassDouble<LaserAbsorptionCoefficient>(ca, "ReferenceTemperature", this, &LaserAbsorptionCoefficient::T0);
  new ClassDouble<LaserAbsorptionCoefficient>(ca, "ReferenceCoefficient", this, &LaserAbsorptionCoefficient::alpha0);
  return ca;
}

//------------------------------------------------------------------------------


LaserData::LaserData() {

  //Rule for laser source: user-specified power time-history overrides (fixed) power, 
  // which overrides (fixed) intensity
  source_intensity = -1.0; //a negative number means not specified
  source_power = -1.0; //a negative number means not specified
  source_power_timehistory_file = "";

  source_distribution = CONSTANT;
  source_radius = 0.0;
  source_beam_waist=0.0;
  source_center_x = source_center_y = source_center_z = 0.0;
  source_dir_x = source_dir_y = source_dir_z = 0.0;
  focusing_angle_degrees = 0.0;
  range = DBL_MAX;

  lmin = 1.0e-12;

  parallel = BALANCED;
  min_cells_per_core = 100;

  source_depth = 0.0;
  alpha = 1.0;
  convergence_tol = 1.0e-5;
  max_iter = 400;
  relax_coeff = 1.0;

}

//------------------------------------------------------------------------------


void LaserData::setup(const char *name, ClassAssigner *father) {

  ClassAssigner *ca = new ClassAssigner(name, 23, father); 

  //Physical Parameters
  new ClassDouble<LaserData>(ca, "SourceIntensity", this, &LaserData::source_intensity);
  new ClassToken<LaserData> (ca, "SourceDistribution", this,
        reinterpret_cast<int LaserData::*>(&LaserData::source_distribution), 2, "Uniform", 0, "Gaussian", 1);
  new ClassDouble<LaserData>(ca, "SourcePower", this, &LaserData::source_power);
  new ClassStr<LaserData>(ca, "SourcePowerTimeHistory", this, &LaserData::source_power_timehistory_file);
  new ClassDouble<LaserData>(ca, "SourceCenterX", this, &LaserData::source_center_x);
  new ClassDouble<LaserData>(ca, "SourceCenterY", this, &LaserData::source_center_y);
  new ClassDouble<LaserData>(ca, "SourceCenterZ", this, &LaserData::source_center_z);
  new ClassDouble<LaserData>(ca, "SourceRadius", this, &LaserData::source_radius);
  new ClassDouble<LaserData>(ca, "SourceBeamWaist", this, &LaserData::source_beam_waist);
  new ClassDouble<LaserData>(ca, "DirectionX", this, &LaserData::source_dir_x);
  new ClassDouble<LaserData>(ca, "DirectionY", this, &LaserData::source_dir_y);
  new ClassDouble<LaserData>(ca, "DirectionZ", this, &LaserData::source_dir_z);
  new ClassDouble<LaserData>(ca, "FocusingAngle", this, &LaserData::focusing_angle_degrees);
  new ClassDouble<LaserData>(ca, "Range", this, &LaserData::range);
  new ClassDouble<LaserData>(ca, "RadianceCutOff", this, &LaserData::lmin);

  abs.setup("AbsorptionCoefficient", ca);

  //Parallelization approach
  new ClassToken<LaserData> (ca, "Parallelization", this,
        reinterpret_cast<int LaserData::*>(&LaserData::parallel), 2, "Original", 0, "Balanced", 1);
  new ClassInt<LaserData>(ca, "NumberOfCellsPerCore", this, &LaserData::min_cells_per_core);

  //Numerical Parameters
  new ClassDouble<LaserData>(ca, "SourceDepth", this, &LaserData::source_depth);
  new ClassDouble<LaserData>(ca, "Alpha", this, &LaserData::alpha);
  new ClassDouble<LaserData>(ca, "ConvergenceTolerance", this, &LaserData::convergence_tol);
  new ClassDouble<LaserData>(ca, "MaxIts", this, &LaserData::max_iter);
  new ClassDouble<LaserData>(ca, "RelaxationCoefficient", this, &LaserData::relax_coeff);

}

//------------------------------------------------------------------------------

AtomicIonizationModel::AtomicIonizationModel()
{
  molar_fraction = -1.0;
  molar_mass = -1.0;
  atomic_number = -1;
  max_charge = 10;
  ionization_energy_filename = "";
  excitation_energy_files_prefix = "";
  excitation_energy_files_suffix = "";
  degeneracy_files_prefix = ""; 
  degeneracy_files_suffix = ""; 
}

//------------------------------------------------------------------------------

Assigner* AtomicIonizationModel::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 9, nullAssigner);

  new ClassDouble<AtomicIonizationModel>(ca, "MolarFraction", this, 
          &AtomicIonizationModel::molar_fraction);
  
  new ClassDouble<AtomicIonizationModel>(ca, "MolarMass", this, 
          &AtomicIonizationModel::molar_mass);
  
  new ClassInt<AtomicIonizationModel>(ca, "AtomicNumber", this, 
          &AtomicIonizationModel::atomic_number);
  
  new ClassInt<AtomicIonizationModel>(ca, "MaxChargeNumber", this, 
          &AtomicIonizationModel::max_charge);
  
  new ClassStr<AtomicIonizationModel>(ca, "IonizationEnergyFile", this, 
          &AtomicIonizationModel::ionization_energy_filename);

  new ClassStr<AtomicIonizationModel>(ca, "ExcitationEnergyFilesPrefix", this, 
          &AtomicIonizationModel::excitation_energy_files_prefix);

  new ClassStr<AtomicIonizationModel>(ca, "ExcitationEnergyFilesSuffix", this, 
          &AtomicIonizationModel::excitation_energy_files_suffix);

  new ClassStr<AtomicIonizationModel>(ca, "DegeneracyFilesPrefix", this, 
          &AtomicIonizationModel::degeneracy_files_prefix);

  new ClassStr<AtomicIonizationModel>(ca, "DegeneracyFilesSuffix", this, 
          &AtomicIonizationModel::degeneracy_files_suffix);

  return ca;

}

//------------------------------------------------------------------------------

MaterialIonizationModel::MaterialIonizationModel()
{
  type = NONE;
  depression = GRIEM; //relevant only to non-ideal Saha
  depression_max = 1.0; //relevant only to non-ideal Saha

  maxIts = 200;
  convergence_tol = 1.0e-5;

  ionization_Tmin = 400; //Kelvin

  partition_evaluation = CUBIC_SPLINE_INTERPOLATION;

  sample_size = 6000;
  sample_Tmin = 10.0; //Kelvin
  sample_Tmax = 1.0e10; //Kelvin
}

//------------------------------------------------------------------------------

Assigner* MaterialIonizationModel::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 11, nullAssigner);

  new ClassToken<MaterialIonizationModel> (ca, "Type", this,
        reinterpret_cast<int MaterialIonizationModel::*>(&MaterialIonizationModel::type), 
        3, "None", 0, "IdealSahaEquation", 1, "NonIdealSahaEquation", 2);

  new ClassToken<MaterialIonizationModel> (ca, "DepressionModel", this,
        reinterpret_cast<int MaterialIonizationModel::*>(&MaterialIonizationModel::depression), 
        4, "None", 0, "Griem", 1, "Ebeling", 2, "GriemFletcher", 3);

  new ClassDouble<MaterialIonizationModel>(ca, "MaxDepression", this, 
        &MaterialIonizationModel::depression_max);

  new ClassToken<MaterialIonizationModel> (ca, "PartitionFunctionEvaluation", this,
        reinterpret_cast<int MaterialIonizationModel::*>(&MaterialIonizationModel::partition_evaluation), 
        3, "OnTheFly", 0, "CubicSplineInterpolation", 1, "LinearInterpolation", 2);

  new ClassDouble<MaterialIonizationModel>(ca, "Tmin", this, 
        &MaterialIonizationModel::ionization_Tmin);

  new ClassDouble<MaterialIonizationModel>(ca, "SampleTmin", this, 
        &MaterialIonizationModel::sample_Tmin);

  new ClassDouble<MaterialIonizationModel>(ca, "SampleTmax", this, 
        &MaterialIonizationModel::sample_Tmax);

  new ClassInt<MaterialIonizationModel>(ca, "SampleSize", this, &MaterialIonizationModel::sample_size);

  new ClassInt<MaterialIonizationModel>(ca, "MaxIts", this, &MaterialIonizationModel::maxIts);

  new ClassDouble<MaterialIonizationModel>(ca, "ConvergenceTolerance", this, 
        &MaterialIonizationModel::convergence_tol);

  elementMap.setup("Element", ca);

  return ca;
}

//------------------------------------------------------------------------------

IonizationData::IonizationData()
{
  planck_constant = 6.62607015e-25;  //unit: (mm^2).g/s  (dim: [energy]/[freq])
  electron_charge = 1.602176634e-19;  //unit: A.s (Coulombs)
  electron_mass = 9.1093837015e-28; //unit: g
  boltzmann_constant = 1.38064852e-14; //unit: (mm^2).g/(s^2*K)  (dim: [energy]/[temperature])
  vacuum_permittivity = 8.8541878128e-24; //unit: (s^4.A^2)/(g*mm^3)  (dim: [electrical capacitance]/[length])
}

//------------------------------------------------------------------------------

void IonizationData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassDouble<IonizationData>(ca, "PlanckConstant", this, &IonizationData::planck_constant);
  new ClassDouble<IonizationData>(ca, "ElectronCharge", this, &IonizationData::electron_charge);
  new ClassDouble<IonizationData>(ca, "ElectronMass", this, &IonizationData::electron_mass);
  new ClassDouble<IonizationData>(ca, "BoltzmannConstant", this, &IonizationData::boltzmann_constant);
  new ClassDouble<IonizationData>(ca, "VacuumPermittivity", this, &IonizationData::vacuum_permittivity);

  materialMap.setup("Material", ca);

}

//------------------------------------------------------------------------------

MaterialVolumes::MaterialVolumes()
{
  filename = "";
  frequency = 10;
  frequency_dt = -1.0;
}

//------------------------------------------------------------------------------

void MaterialVolumes::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassStr<MaterialVolumes>(ca, "FileName", this, &MaterialVolumes::filename);
  new ClassInt<MaterialVolumes>(ca, "Frequency", this, &MaterialVolumes::frequency);
  new ClassDouble<MaterialVolumes>(ca, "TimeInterval", this, &MaterialVolumes::frequency_dt);

}

//------------------------------------------------------------------------------

OutputData::OutputData()
{
  prefix = "";
  solution_filename_base = "solution";

  frequency = 0;
  frequency_dt = -1.0;

  density = OFF;
  velocity = OFF;
  pressure = OFF;
  materialid = OFF;
  temperature = OFF;
  delta_temperature = OFF;
  internal_energy = OFF;
  delta_internal_energy = OFF;
  laser_radiance = OFF;
  reference_map = OFF;
  principal_elastic_stresses = OFF;

  levelset0 = OFF;
  levelset1 = OFF;
  levelset2 = OFF;
  levelset3 = OFF;
  levelset4 = OFF;

  for(int i=0; i<MAXLS; i++)
    levelset[i] = OFF;

  mean_charge = OFF;
  heavy_particles_density = OFF;
  electron_density = OFF;
  max_charge_number = 3;
  molar_fractions0 = OFF; 
  molar_fractions1 = OFF; 
  molar_fractions2 = OFF; 
  molar_fractions3 = OFF; 
  molar_fractions4 = OFF; 

  for(int i=0; i<MAXSPECIES; i++)
    molar_fractions[i] = OFF;

  mesh_filename = "";

  mesh_partition = "";

  verbose = LOW;
}

//------------------------------------------------------------------------------

void OutputData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 26+MAXLS+MAXSPECIES, father);

  new ClassStr<OutputData>(ca, "Prefix", this, &OutputData::prefix);
  new ClassStr<OutputData>(ca, "Solution", this, &OutputData::solution_filename_base);

  new ClassInt<OutputData>(ca, "Frequency", this, &OutputData::frequency);
  new ClassDouble<OutputData>(ca, "TimeInterval", this, &OutputData::frequency_dt);

  new ClassToken<OutputData>(ca, "Density", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::density), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "Velocity", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::velocity), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "Pressure", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::pressure), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "MaterialID", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::materialid), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "Temperature", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::temperature), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "DeltaTemperature", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::delta_temperature), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "InternalEnergyPerUnitMass", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::internal_energy), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "DeltaInternalEnergyPerUnitMass", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::delta_internal_energy), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "LaserRadiance", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::laser_radiance), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "ReferenceMap", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::reference_map), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "PrincipalElasticStresses", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::principal_elastic_stresses), 2,
                             "Off", 0, "On", 1);

  new ClassToken<OutputData>(ca, "LevelSet0", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::levelset0), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "LevelSet1", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::levelset1), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "LevelSet2", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::levelset2), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "LevelSet3", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::levelset3), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "LevelSet4", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::levelset4), 2,
                             "Off", 0, "On", 1);

  new ClassToken<OutputData>(ca, "MeanCharge", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::mean_charge), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "HeavyParticlesDensity", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::heavy_particles_density), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "ElectronDensity", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::electron_density), 2,
                             "Off", 0, "On", 1);

  new ClassInt<OutputData>(ca, "MaxChargeNumber", this, &OutputData::max_charge_number);

  new ClassToken<OutputData>(ca, "MolarFractionsElement0", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::molar_fractions0), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "MolarFractionsElement1", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::molar_fractions1), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "MolarFractionsElement2", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::molar_fractions2), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "MolarFractionsElement3", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::molar_fractions3), 2,
                             "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "MolarFractionsElement4", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::molar_fractions4), 2,
                             "Off", 0, "On", 1);

  new ClassStr<OutputData>(ca, "MeshInformation", this, &OutputData::mesh_filename);

  new ClassStr<OutputData>(ca, "MeshPartition", this, &OutputData::mesh_partition);

  new ClassToken<OutputData>(ca, "VerboseScreenOutput", this,
                             reinterpret_cast<int OutputData::*>(&OutputData::verbose), 3,
                             "Low", 0, "Medium", 1, "High", 2);

  probes.setup("Probes", ca);

  linePlots.setup("LinePlot", ca);

  planePlots.setup("CutPlane", ca);

  materialVolumes.setup("MaterialVolumes", ca);

}

//------------------------------------------------------------------------------

ProbeNode::ProbeNode() {
  locationX = locationY = locationZ = -1.0e20;
}

//------------------------------------------------------------------------------

Assigner* ProbeNode::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 3, nullAssigner);

  new ClassDouble<ProbeNode>(ca, "X",this,&ProbeNode::locationX);
  new ClassDouble<ProbeNode>(ca, "Y",this,&ProbeNode::locationY);
  new ClassDouble<ProbeNode>(ca, "Z",this,&ProbeNode::locationZ);

  return ca;
}

//------------------------------------------------------------------------------

Probes::Probes() {

  frequency = -100;
  frequency_dt = -1;

  density = "";
  pressure = "";
  temperature = "";
  delta_temperature = "";
  velocity_x = "";
  velocity_y = "";
  velocity_z = "";
  materialid = "";
  laser_radiance = "";
  levelset0 = "";
  levelset1 = "";
  levelset2 = "";
  levelset3 = "";
  levelset4 = "";
  ionization_result = "";
  reference_map = "";
  principal_elastic_stresses = "";

}

//------------------------------------------------------------------------------

void Probes::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 20, father);

  new ClassInt<Probes>(ca, "Frequency", this, &Probes::frequency);
  new ClassDouble<Probes>(ca, "TimeInterval", this, &Probes::frequency_dt);
  new ClassStr<Probes>(ca, "Density", this, &Probes::density);
  new ClassStr<Probes>(ca, "Pressure", this, &Probes::pressure);
  new ClassStr<Probes>(ca, "Temperature", this, &Probes::temperature);
  new ClassStr<Probes>(ca, "DeltaTemperature", this, &Probes::delta_temperature);
  new ClassStr<Probes>(ca, "VelocityX", this, &Probes::velocity_x);
  new ClassStr<Probes>(ca, "VelocityY", this, &Probes::velocity_y);
  new ClassStr<Probes>(ca, "VelocityZ", this, &Probes::velocity_z);
  new ClassStr<Probes>(ca, "MaterialID", this, &Probes::materialid);
  new ClassStr<Probes>(ca, "LaserRadiance", this, &Probes::laser_radiance);
  new ClassStr<Probes>(ca, "LevelSet0", this, &Probes::levelset0);
  new ClassStr<Probes>(ca, "LevelSet1", this, &Probes::levelset1);
  new ClassStr<Probes>(ca, "LevelSet2", this, &Probes::levelset2);
  new ClassStr<Probes>(ca, "LevelSet3", this, &Probes::levelset3);
  new ClassStr<Probes>(ca, "LevelSet4", this, &Probes::levelset4);
  new ClassStr<Probes>(ca, "IonizationResult", this, &Probes::ionization_result);
  new ClassStr<Probes>(ca, "ReferenceMap", this, &Probes::reference_map);
  new ClassStr<Probes>(ca, "PrincipalElasticStresses", this, &Probes::principal_elastic_stresses);

  myNodes.setup("Node", ca);

}

//------------------------------------------------------------------------------

LinePlot::LinePlot() {

  x0 = y0 = z0 = 0.0;
  x1 = y1 = z1 = 0.0;

  numPoints = 0;
  frequency = -100;
  frequency_dt = -1.0;

  filename_base = "";
}

//------------------------------------------------------------------------------

Assigner* LinePlot::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 10, nullAssigner);

  new ClassStr<LinePlot>(ca, "FileName", this, &LinePlot::filename_base);

  new ClassInt<LinePlot>(ca, "NumberOfPoints", this, &LinePlot::numPoints);
  new ClassInt<LinePlot>(ca, "Frequency", this, &LinePlot::frequency);
  new ClassDouble<LinePlot>(ca, "TimeInterval", this, &LinePlot::frequency_dt);

  new ClassDouble<LinePlot>(ca, "X0", this, &LinePlot::x0);
  new ClassDouble<LinePlot>(ca, "Y0", this, &LinePlot::y0);
  new ClassDouble<LinePlot>(ca, "Z0", this, &LinePlot::z0);
  new ClassDouble<LinePlot>(ca, "Xmax", this, &LinePlot::x1);
  new ClassDouble<LinePlot>(ca, "Ymax", this, &LinePlot::y1);
  new ClassDouble<LinePlot>(ca, "Zmax", this, &LinePlot::z1);

  return ca;

}

//------------------------------------------------------------------------------

PlanePlot::PlanePlot() {

  x0 = y0 = z0 = 0.0;
  normal_x = 0.0;
  normal_y = 0.0;
  normal_z = 1.0;

  frequency = -100;
  frequency_dt = -1.0;

  mesh = "";

  density = "";
  pressure = "";
  temperature = "";
  delta_temperature = "";
  velocity = "";
  materialid = "";
  laser_radiance = "";
  levelset0 = "";
  levelset1 = "";
  levelset2 = "";
  levelset3 = "";
  levelset4 = "";
  ionization_result = "";

}

//------------------------------------------------------------------------------

Assigner* PlanePlot::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 22, nullAssigner);

  new ClassInt<PlanePlot>(ca, "Frequency", this, &PlanePlot::frequency);
  new ClassDouble<PlanePlot>(ca, "TimeInterval", this, &PlanePlot::frequency_dt);

  new ClassDouble<PlanePlot>(ca, "X0", this, &PlanePlot::x0);
  new ClassDouble<PlanePlot>(ca, "Y0", this, &PlanePlot::y0);
  new ClassDouble<PlanePlot>(ca, "Z0", this, &PlanePlot::z0);
  new ClassDouble<PlanePlot>(ca, "Normal_x", this, &PlanePlot::normal_x);
  new ClassDouble<PlanePlot>(ca, "Normal_y", this, &PlanePlot::normal_y);
  new ClassDouble<PlanePlot>(ca, "Normal_z", this, &PlanePlot::normal_z);

  new ClassStr<PlanePlot>(ca, "Mesh", this, &PlanePlot::mesh);

  new ClassStr<PlanePlot>(ca, "Density", this, &PlanePlot::density);
  new ClassStr<PlanePlot>(ca, "Pressure", this, &PlanePlot::pressure);
  new ClassStr<PlanePlot>(ca, "Temperature", this, &PlanePlot::temperature);
  new ClassStr<PlanePlot>(ca, "DeltaTemperature", this, &PlanePlot::delta_temperature);
  new ClassStr<PlanePlot>(ca, "Velocity", this, &PlanePlot::velocity);
  new ClassStr<PlanePlot>(ca, "MaterialID", this, &PlanePlot::materialid);
  new ClassStr<PlanePlot>(ca, "LaserRadiance", this, &PlanePlot::laser_radiance);
  new ClassStr<PlanePlot>(ca, "LevelSet0", this, &PlanePlot::levelset0);
  new ClassStr<PlanePlot>(ca, "LevelSet1", this, &PlanePlot::levelset1);
  new ClassStr<PlanePlot>(ca, "LevelSet2", this, &PlanePlot::levelset2);
  new ClassStr<PlanePlot>(ca, "LevelSet3", this, &PlanePlot::levelset3);
  new ClassStr<PlanePlot>(ca, "LevelSet4", this, &PlanePlot::levelset4);
  new ClassStr<PlanePlot>(ca, "IonizationResult", this, &PlanePlot::ionization_result);

  return ca;

}

//------------------------------------------------------------------------------


EmbeddedSurfaceData::EmbeddedSurfaceData()
{
  provided_by_another_solver = NO;

  surface_thickness = 1.0e-8;

  // force calculation
  gauss_points_lofting = 0.0;
  internal_pressure = 0.0;
  quadrature = ONE_POINT;
  twoD_to_threeD = RADIAL_BASIS;

  filename = "";
  type = None;
  thermal  = Adiabatic;
  wall_temperature = 300.0; //!< Kelvin
  heat_source = 0.0;
  dynamics_calculator = "";
  force_calculator = "";


  conRec_depth = 0.0;
}

//------------------------------------------------------------------------------

Assigner *EmbeddedSurfaceData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 14, nullAssigner);

  new ClassToken<EmbeddedSurfaceData> (ca, "SurfaceProvidedByAnotherSolver", this,
     reinterpret_cast<int EmbeddedSurfaceData::*>(&EmbeddedSurfaceData::provided_by_another_solver), 2,
     "No", 0, "Yes", 1);

  new ClassDouble<EmbeddedSurfaceData>(ca, "SurfaceThickness", this, 
                                      &EmbeddedSurfaceData::surface_thickness);

  new ClassStr<EmbeddedSurfaceData>(ca, "MeshFile", this, &EmbeddedSurfaceData::filename);


  new ClassToken<EmbeddedSurfaceData> (ca, "GaussQuadrature", this,
     reinterpret_cast<int EmbeddedSurfaceData::*>(&EmbeddedSurfaceData::quadrature), 8,
     "None", 0, "OnePoint", 1, "ThreePoint", 2, "ThreePoints", 2,
     "FourPoint", 3, "FourPoints", 3, "SixPoint", 4, "SixPoints", 4);

  new ClassDouble<EmbeddedSurfaceData>(ca, "GaussPointsLofting", this, &EmbeddedSurfaceData::gauss_points_lofting);

  new ClassDouble<EmbeddedSurfaceData>(ca, "InternalPressure", this, &EmbeddedSurfaceData::internal_pressure);

  new ClassToken<EmbeddedSurfaceData> (ca, "TwoDimensionalToThreeDimensionalMapping", this,
      reinterpret_cast<int EmbeddedSurfaceData::*>(&EmbeddedSurfaceData::twoD_to_threeD), 2,
      "RadialBasisInterpolation", 0, "NearestNeighbor", 1);

  new ClassToken<EmbeddedSurfaceData> (ca, "BoundaryCondition", this,
     reinterpret_cast<int EmbeddedSurfaceData::*>(&EmbeddedSurfaceData::type), 6,
     "None", 0, "Wall", 1, "Symmetry", 2, "DirectState", 3, "MassFlow", 4, "PorousWall", 5);

  new ClassToken<EmbeddedSurfaceData> (ca, "ThermalBoundaryCondition", this,
     reinterpret_cast<int EmbeddedSurfaceData::*>(&EmbeddedSurfaceData::thermal), 3,
     "Adiabatic", 0, "Isothermal", 1, "Source", 2);

  new ClassDouble<EmbeddedSurfaceData>(ca, "Temperature", this, &EmbeddedSurfaceData::wall_temperature);

  new ClassDouble<EmbeddedSurfaceData>(ca, "HeatSource", this, &EmbeddedSurfaceData::heat_source);

  new ClassStr<EmbeddedSurfaceData>(ca, "UserDefinedDynamicsCalculator", this, 
                                    &EmbeddedSurfaceData::dynamics_calculator);

  new ClassStr<EmbeddedSurfaceData>(ca, "UserDefinedForceCalculator", this, 
                                    &EmbeddedSurfaceData::force_calculator);

  new ClassDouble<EmbeddedSurfaceData>(ca, "ConstantReconstructionDepth", 
                                       this, &EmbeddedSurfaceData::conRec_depth);

  output.setup("Output", ca); //there is another "Output", must provide "ca" to distinguish

  return ca;
}

//------------------------------------------------------------------------------

EmbeddedSurfacesData::EmbeddedSurfacesData()
{

}

//------------------------------------------------------------------------------

void EmbeddedSurfacesData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  surfaces.setup("Surface", ca);
}

//------------------------------------------------------------------------------

EmbeddedBoundaryMethodData::EmbeddedBoundaryMethodData()
{
  riemann_normal = MESH;
  recon = CONSTANT;
}

//------------------------------------------------------------------------------

void EmbeddedBoundaryMethodData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  embed_surfaces.setup("EmbeddedSurfaces");

  new ClassToken<EmbeddedBoundaryMethodData>
    (ca, "RiemannNormal", this,
     reinterpret_cast<int EmbeddedBoundaryMethodData::*>(&EmbeddedBoundaryMethodData::riemann_normal), 
     3, "EmbeddedSurface", 0, "Mesh", 1, "Average", 2);

  new ClassToken<EmbeddedBoundaryMethodData>
    (ca, "ReconstructionAtInterface", this,
     reinterpret_cast<int EmbeddedBoundaryMethodData::*>(&EmbeddedBoundaryMethodData::recon),
     2, "Constant", 0, "Linear", 1);
}

//------------------------------------------------------------------------------

AerosCouplingData::AerosCouplingData()
{
  fsi_algo = NONE;
}

//------------------------------------------------------------------------------

void AerosCouplingData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<AerosCouplingData> (ca, "FSIAlgorithm", this,
     reinterpret_cast<int AerosCouplingData::*>(&AerosCouplingData::fsi_algo), 4,
     "None", 0, "ByAeroS", 1, "C0", 2, "A6", 3);
}

//------------------------------------------------------------------------------

AerofCouplingData::AerofCouplingData()
{
  type = NONE;
}

//------------------------------------------------------------------------------

void AerofCouplingData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<AerofCouplingData> (ca, "Type", this,
     reinterpret_cast<int AerofCouplingData::*>(&AerofCouplingData::type), 2,
     "None", 0, "OversetGrids", 1);
}

//------------------------------------------------------------------------------

M2CTwinningData::M2CTwinningData()
{
  type = NONE;
}

//------------------------------------------------------------------------------

void M2CTwinningData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<M2CTwinningData> (ca, "Type", this,
     reinterpret_cast<int M2CTwinningData::*>(&M2CTwinningData::type), 2,
     "None", 0, "OversetGrids", 1);
}

//------------------------------------------------------------------------------


ConcurrentProgramsData::ConcurrentProgramsData()
{

}

//------------------------------------------------------------------------------

void ConcurrentProgramsData::setup(const char *name, ClassAssigner *father)
{
  //ClassAssigner *ca = new ClassAssigner(name, 3, father);
  new ClassAssigner(name, 3, father);
  aeros.setup("AeroS");
  aerof.setup("AeroF");
  m2c_twin.setup("M2CTwin");
} 

//------------------------------------------------------------------------------

LagrangianMeshOutputData::LagrangianMeshOutputData()
{
  frequency = 0;
  frequency_dt = -1.0;

  prefix = "";

  orig_config = "";
  disp = "";
  sol  = "";

  wetting_output_filename = "";
}

//------------------------------------------------------------------------------

void LagrangianMeshOutputData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 7, father);
  
  new ClassInt<LagrangianMeshOutputData>(ca, "Frequency", this, &LagrangianMeshOutputData::frequency);
  new ClassDouble<LagrangianMeshOutputData>(ca, "TimeInterval", this, &LagrangianMeshOutputData::frequency_dt);

  new ClassStr<LagrangianMeshOutputData>(ca, "Prefix", this, &LagrangianMeshOutputData::prefix);
  new ClassStr<LagrangianMeshOutputData>(ca, "Mesh", this, &LagrangianMeshOutputData::orig_config);

  new ClassStr<LagrangianMeshOutputData>(ca, "Displacement", this, &LagrangianMeshOutputData::disp);
  new ClassStr<LagrangianMeshOutputData>(ca, "Solution", this, &LagrangianMeshOutputData::sol);

  new ClassStr<LagrangianMeshOutputData>(ca, "ContactSurfaceOutput", this,
                                         &LagrangianMeshOutputData::wetting_output_filename);

}

//------------------------------------------------------------------------------

TransientInputData::TransientInputData()
{
  metafile = "";
  snapshot_file_prefix = "";
  snapshot_file_suffix = "";

  basis = INVERSE_MULTIQUADRIC;
  numPoints = 8;
}

//------------------------------------------------------------------------------

void TransientInputData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 6, father);
  
  new ClassStr<TransientInputData>(ca, "MetaFile", this, &TransientInputData::metafile);

  new ClassStr<TransientInputData>(ca, "SnapshotFilePrefix", this, 
          &TransientInputData::snapshot_file_prefix);

  new ClassStr<TransientInputData>(ca, "SnapshotFileSuffix", this, 
          &TransientInputData::snapshot_file_suffix);

  new ClassToken<TransientInputData> (ca, "SpatialInterpolationBasis", this,
     reinterpret_cast<int TransientInputData::*>(&TransientInputData::basis), 4,
     "Multiquadric", 0, "InverseMultiquadric", 1, "ThinPlateSpline", 2, "Gaussian", 3);

  new ClassInt<TransientInputData>(ca, "NumberOfBasisPoints", this, &TransientInputData::numPoints);

  output.setup("Output", ca); //there is another "Output", must provide "ca" to distinguish
}

//------------------------------------------------------------------------------

TerminalVisualizationData::TerminalVisualizationData()
{
  colormap = TURBO;
  plane = NONE;
  filename = "";
  variable = PRESSURE;
  coordinate     = DBL_MAX;
  horizontal_min = -DBL_MAX;
  horizontal_max = DBL_MAX;
  vertical_min   = -DBL_MAX;
  vertical_max   = DBL_MAX;
  dx = DBL_MAX;

  frequency = 0;
  frequency_dt = -1.0;
  frequency_clocktime = 60; //seconds
  pause = 2.0; //seconds 
}

//------------------------------------------------------------------------------

void TerminalVisualizationData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 14, father);

  new ClassToken<TerminalVisualizationData> (ca, "ColorMap", this,
      reinterpret_cast<int TerminalVisualizationData::*>(&TerminalVisualizationData::colormap), 2,
      "GrayScale", 0, "Turbo", 1);

  new ClassToken<TerminalVisualizationData> (ca, "Plane", this,
      reinterpret_cast<int TerminalVisualizationData::*>(&TerminalVisualizationData::plane), 4,
      "None", 0, "YZ", 1, "XZ", 2, "XY", 3);

  new ClassDouble<TerminalVisualizationData>(ca, "Coordinate", this, 
      &TerminalVisualizationData::coordinate);

  new ClassStr<TerminalVisualizationData>(ca, "FileName", this, 
      &TerminalVisualizationData::filename);

  new ClassToken<TerminalVisualizationData> (ca, "Variable", this,
      reinterpret_cast<int TerminalVisualizationData::*>(&TerminalVisualizationData::variable), 9,
      "Density", 0, "Velocity", 1, "Pressure", 2, "Temperature", 3, "MaterialID", 4, "LaserRadiance", 5,
      "LevelSet0", 6, "LevelSet1", 7, "MeanCharge", 8);

  new ClassDouble<TerminalVisualizationData>(ca, "HorizontalMin", this, 
      &TerminalVisualizationData::horizontal_min);

  new ClassDouble<TerminalVisualizationData>(ca, "HorizontalMax", this, 
      &TerminalVisualizationData::horizontal_max);

  new ClassDouble<TerminalVisualizationData>(ca, "VerticalMin", this, 
      &TerminalVisualizationData::vertical_min);

  new ClassDouble<TerminalVisualizationData>(ca, "VerticalMax", this, 
      &TerminalVisualizationData::vertical_max);

  new ClassDouble<TerminalVisualizationData>(ca, "Resolution", this, 
      &TerminalVisualizationData::dx);

  new ClassInt<TerminalVisualizationData>(ca, "Frequency", this, 
      &TerminalVisualizationData::frequency);

  new ClassDouble<TerminalVisualizationData>(ca, "TimeInterval", this, 
      &TerminalVisualizationData::frequency_dt);

  new ClassDouble<TerminalVisualizationData>(ca, "ClockTimeInterval", this, 
      &TerminalVisualizationData::frequency_clocktime);

  new ClassDouble<TerminalVisualizationData>(ca, "Pause", this, 
      &TerminalVisualizationData::pause);

}

//------------------------------------------------------------------------------

ReferenceMapData::ReferenceMapData()
{
  fd = UPWIND_CENTRAL_3;
}

//------------------------------------------------------------------------------

void ReferenceMapData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<ReferenceMapData> (ca, "FiniteDifference", this,
     reinterpret_cast<int ReferenceMapData::*>(&ReferenceMapData::fd), 2,
     "None", 0, "ThirdOrderUpwind", 1);
}

//------------------------------------------------------------------------------

EOSTabulationData::EOSTabulationData()
{
  materialid = -1;
  filename = ""; 
  output = PRESSURE;
  xvar = DENSITY;
  yvar = SPECIFIC_INTERNAL_ENERGY;
  x0 = xmax = y0 = ymax = 1.0;
  Nx = Ny = 100;
}

//------------------------------------------------------------------------------

Assigner *EOSTabulationData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 11, nullAssigner);

  new ClassInt<EOSTabulationData>(ca, "MaterialID", this, &EOSTabulationData::materialid);

  new ClassStr<EOSTabulationData>(ca, "OutputFile", this, &EOSTabulationData::filename);

  new ClassToken<EOSTabulationData> (ca, "TabulatedVariable", this,
          reinterpret_cast<int EOSTabulationData::*>(&EOSTabulationData::output), 10,
          "Pressure", 0, 
          "SpecificInternalEnergy", 1, "InternalEnergyPerUnitMass", 1,
          "Density", 2, "PressureDerivativeEnergy", 3, 
          "GruneisenParameter", 4, "PressureDerivativeDensity", 5,
          "BulkModulus", 6, "Temperature", 7, 
          "SpecificEnthalpy", 8, "EnthalpyPerUnitMass", 8);

  new ClassToken<EOSTabulationData> (ca, "VariableX", this,
          reinterpret_cast<int EOSTabulationData::*>(&EOSTabulationData::xvar), 5,
          "Density", 2, 
          "SpecificInternalEnergy", 1, "InternalEnergyPerUnitMass", 1,
          "Pressure", 0, "Temperature", 7);

  new ClassToken<EOSTabulationData> (ca, "VariableY", this,
          reinterpret_cast<int EOSTabulationData::*>(&EOSTabulationData::yvar), 5,
          "Density", 2, 
          "SpecificInternalEnergy", 1, "InternalEnergyPerUnitMass", 1,
          "Pressure", 0, "Temperature", 7);


  new ClassDouble<EOSTabulationData>(ca, "X0", this, &EOSTabulationData::x0);
  new ClassDouble<EOSTabulationData>(ca, "Xmax", this, &EOSTabulationData::xmax);
  new ClassDouble<EOSTabulationData>(ca, "Y0", this, &EOSTabulationData::y0);
  new ClassDouble<EOSTabulationData>(ca, "Ymax", this, &EOSTabulationData::ymax);

  new ClassInt<EOSTabulationData>(ca, "NumberOfPointsX", this, &EOSTabulationData::Nx);
  new ClassInt<EOSTabulationData>(ca, "NumberOfPointsY", this, &EOSTabulationData::Ny);

  return ca;
}

//------------------------------------------------------------------------------

SpecialToolsData::SpecialToolsData()
{
  type = NONE;
}

//------------------------------------------------------------------------------

void SpecialToolsData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassToken<SpecialToolsData> (ca, "Type", this,
     reinterpret_cast<int SpecialToolsData::*>(&SpecialToolsData::type), 3,
     "None", 0, "DynamicLoadCalculation", 1, "EquationOfStateTabulation", 2);

  transient_input.setup("TransientInputData");
  eos_tabulationMap.setup("EquationOfStateTable", ca);
} 

//------------------------------------------------------------------------------

IoData::IoData(int argc, char** argv)
{
  //Should NOT call functions in Utils (e.g., print(), exit_mpi()) because the
  //M2C communicator may have not been properly set up.
  readCmdLine(argc, argv);
  readCmdFile();
}

//------------------------------------------------------------------------------

void IoData::readCmdLine(int argc, char** argv)
{
  if(argc==1) {
    fprintf(stdout,"\033[0;31m*** Error: Input file not provided!\n\033[0m");
    exit(-1);
  }
  cmdFileName = argv[1];
}

//------------------------------------------------------------------------------

void IoData::readCmdFile()
{
  extern FILE *yyCmdfin;
  extern int yyCmdfparse();

  setupCmdFileVariables();
//  cmdFilePtr = freopen(cmdFileName, "r", stdin);
  yyCmdfin = cmdFilePtr = fopen(cmdFileName, "r");

  if (!cmdFilePtr) {
    fprintf(stdout,"\033[0;31m*** Error: could not open \'%s\'\n\033[0m", cmdFileName);
    exit(-1);
  }

  int error = yyCmdfparse();
  if (error) {
    fprintf(stdout,"\033[0;31m*** Error: command file contained parsing errors.\n\033[0m");
    exit(error);
  }
  fclose(cmdFilePtr);
}

//------------------------------------------------------------------------------
// This function is supposed to be called after creating M2C communicator. So, 
// functions in Utils can be used.
void IoData::finalize()
{
  //Check spatial domain (for spherical and cylindrical)
  mesh.check();

  //READ ADDITIONAL FILES
  if(strcmp(ic.user_specified_ic, ""))
    ic.readUserSpecifiedIC(); 

  //assign default initial state to farfield/inlet b.c. if user did not specify a
  //default (TODO: This logic would fail if user-specified default i.c. is exactly the same
  //as the default StateVariable(). But this is highly unlikely...)
  if(ic.default_ic == StateVariable() &&
     !(bc.inlet == ic.default_ic))
    ic.default_ic = bc.inlet;

  //Set dummy_state (except material id)
  if(fabs(eqs.dummy_state.density-1.0e-6)<1e-12 &&
     eqs.dummy_state.velocity_x == 0.0 && eqs.dummy_state.velocity_y == 0.0 &&
     eqs.dummy_state.velocity_z == 0.0 && eqs.dummy_state.pressure   == 0.0 &&
     eqs.dummy_state.temperature == 0.0 && 
     eqs.dummy_state.internal_energy_per_mass == 0.0) {//user did not specify dummy_state
    eqs.dummy_state = ic.default_ic;
    eqs.dummy_state.materialid = 0; //Not dummy_state's materialid. Shouldn't be used
  }

  //FIX Levelset output (TODO: need a better way...)
  output.levelset[0] = output.levelset0;
  output.levelset[1] = output.levelset1;
  output.levelset[2] = output.levelset2;
  output.levelset[3] = output.levelset3;
  output.levelset[4] = output.levelset4;

  //FIX elemental molar fractions (TODO: need a better way...)
  output.molar_fractions[0] = output.molar_fractions0;
  output.molar_fractions[1] = output.molar_fractions1;
  output.molar_fractions[2] = output.molar_fractions2;
  output.molar_fractions[3] = output.molar_fractions3;
  output.molar_fractions[4] = output.molar_fractions4;
}

//------------------------------------------------------------------------------

void IoData::setupCmdFileVariables()
{

  concurrent.setup("ConcurrentPrograms");

  ebm.setup("EmbeddedBoundary");
  ebm.setup("EmbeddedBoundaryMethod");

  eqs.setup("Equations");
  eqs.setup("NavierStokesEquations");

  ic.setup("InitialCondition");

  bc.setup("BoundaryConditions");

  mesh.setup("Mesh");

  schemes.setup("Space");

  exact_riemann.setup("ExactRiemannSolution");

  laser.setup("Laser");

  ion.setup("Ionization");

  ts.setup("Time");

  multiphase.setup("MultiPhase");

  refmap.setup("ReferenceMap");

  output.setup("Output");

  special_tools.setup("SpecialTools");

  terminal_visualization.setup("TerminalVisualization");

}

//------------------------------------------------------------------------------
