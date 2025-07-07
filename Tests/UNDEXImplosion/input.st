// Wentao Ma et al., IJSS (2022), Mesh Pair #5
// Units: g, mm, s

under ConcurrentPrograms {
  under AeroS {
    FSIAlgorithm = ByAeroS;
  }
}

under Mesh {
  Type = ThreeDimensional;
  X0   = -0.2;
  Xmax =  0.2;
  Y0   = -600.0;
  Ymax =  600.0;
  Z0   = -600.0;
  Zmax =  600.0;

  NumberOfCellsX = 1;

  under ControlPointY[0] {Coordinate =-600.0; CellWidth  = 20.0;}
  under ControlPointY[1] {Coordinate =-160.0; CellWidth  = 3.0;}
  under ControlPointY[2] {Coordinate = -90.0; CellWidth  = 0.7;}
  under ControlPointY[3] {Coordinate = -60.0; CellWidth  = 0.3;}
  under ControlPointY[4] {Coordinate =  30.0; CellWidth  = 0.3;}
  under ControlPointY[5] {Coordinate =  60.0; CellWidth  = 0.7;}
  under ControlPointY[6] {Coordinate = 130.0; CellWidth  = 3.0;}
  under ControlPointY[7] {Coordinate = 600.0; CellWidth  = 20.0;}

  under ControlPointZ[0] {Coordinate =-600.0; CellWidth  = 20.0;}
  under ControlPointZ[1] {Coordinate =-150.0; CellWidth  = 3.0;}
  under ControlPointZ[2] {Coordinate = -70.0; CellWidth  = 0.7;}
  under ControlPointZ[3] {Coordinate = -30.0; CellWidth  = 0.3;}
  under ControlPointZ[4] {Coordinate =  30.0; CellWidth  = 0.3;}
  under ControlPointZ[5] {Coordinate =  70.0; CellWidth  = 0.7;}
  under ControlPointZ[6] {Coordinate = 150.0; CellWidth  = 3.0;}
  under ControlPointZ[7] {Coordinate = 600.0; CellWidth  = 20.0;}

  BoundaryConditionX0   = Symmetry;
  BoundaryConditionXmax = Symmetry;
  BoundaryConditionY0   = Farfield;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Farfield;
  BoundaryConditionZmax = Farfield;
}

under EmbeddedBoundaryMethod {
  RiemannNormal = Mesh; //Average;
  ReconstructionAtInterface = Constant;
  under EmbeddedSurfaces {
    under Surface[0] {
      SurfaceProvidedByAnotherSolver = Yes;
      SurfaceThickness = 1.0e-6;
      GaussQuadrature = OnePoint;
      GaussPointsLofting = 0.2;
      InternalPressure = 0.0;
      ConstantReconstructionDepth = 0.3;
      under Output {
        //Frequency = 1;
        TimeInterval = 1.0e-5;
        Prefix = "results/";
        Mesh = "surface.top";
        Displacement = "surface_disp.txt";
        Solution = "surface_load.txt";
        ContactSurfaceOutput = "surface_contact.top";
      }
    }
  }
}

under Equations {
  under Material[0] { //water
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 7.15;
      PressureConstant = 2.89e8;
    }
  }

  under Material[1] { //air
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
    }
  }

  under Material[2] { //air
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
    }
  }
}

under InitialCondition {
  
  under GeometricEntities {

    under Sphere[0] { //explosion bubble
      Center_x = 0;
      Center_y = -30.0;
      Center_z = 0;
      Radius = 2.5;
      under InitialState {
        MaterialID = 1;
        Density = 5.0e-5;
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 25.0e6;
      }
    }

    under Point[0] { //inside cylinder
      X = 0.0;
      Y = 0.0;
      Z = 0.0;
      under InitialState {
        MaterialID = 2;
        Density = 1.225e-6;
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.0e5;
      }
    }
  }
}

under BoundaryConditions {
  under Farfield {
    MaterialID = 0;
    Density = 0.001;
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e6; //hydrostatic pressure
  }
}

under Space {

  under NavierStokes {
    Flux = LocalLaxFriedrichs; //HLLC; //LocalLaxFriedrichs; //Roe; 
    under Reconstruction {
      Type = Linear;    
      VariableType = ConservativeCharacteristic;
//      SlopeNearInterface = Zero;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.2;
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    Solver = FiniteDifference;
    //Bandwidth = 7;

    BoundaryConditionX0   = ZeroNeumann;
    BoundaryConditionXmax = ZeroNeumann;
    BoundaryConditionY0   = LinearExtrapolation;
    BoundaryConditionYmax = LinearExtrapolation;
    BoundaryConditionZ0   = LinearExtrapolation;
    BoundaryConditionZmax = LinearExtrapolation;

    under Reinitialization {
      Frequency = 8000;
      MaxIts = 400;
    }
  }
}

under MultiPhase {
  Flux = Numerical;
  ReconstructionAtInterface = Constant;
  PhaseChange = RiemannSolution;
  RiemannNormal = LevelSet;
  //LevelSetCorrectionFrequency = 100;
  //ConstantReconstructionDepth = 0.8; //mm
}

under Time {
  Type = Explicit;
  MaxTime = 4.0e-3;
  CFL    = 0.3;
  under Explicit {
    Type = RungeKutta3; 
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  //Frequency = 1;
  TimeInterval = 1.0e-5;
  Density = On;
  Velocity = On;
  Pressure = On;
  LevelSet0 = On;
  MaterialID = On;

  //MeshPartition = "partition";

  under Probes {
    TimeInterval = 2.0e-6;
    Pressure = "pressure_probes.txt";
    VelocityX = "vx_probes.txt";
    VelocityY = "vy_probes.txt";
    VelocityZ = "vz_probes.txt";
    Density  = "density_probes.txt";
    under Node[0] {X = 0.0; Y = -23.3342; Z = 23.3342;}
    under Node[1] {X = 0.0; Y = 0.0; Z = 32.9995;}
    under Node[2] {X = 0.0; Y = 32.9998; Z = 0.0;}
  }

  VerboseScreenOutput = Medium;
}

/*
under TerminalVisualization {
  Plane = YZ;
  Variable = Pressure;
  ClockTimeInterval = 60; //seconds
  HorizontalMin = -40.0;
  HorizontalMax = 20.0;
  VerticalMin = -20.0;
  VerticalMax = 20.0;
}
*/
