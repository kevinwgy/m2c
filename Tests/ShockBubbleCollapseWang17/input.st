//units: g, mm, s, K

under Mesh {

  Type = Cylindrical;
  X0   = -5.0;
  Xmax = 0.05;
  Y0   = 0.0;
  Ymax = 2.0;
  Z0   = -0.001;
  Zmax = 0.001;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = SlipWall;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = SlipWall;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate = -5.0;  CellWidth = 0.1;}
  under ControlPointX[1] {Coordinate = -4.0;  CellWidth = 0.1;}
  under ControlPointX[2] {Coordinate = -1.0;  CellWidth = 0.006;}
  under ControlPointX[3] {Coordinate = -0.3;  CellWidth = 0.0015;}
  under ControlPointX[4] {Coordinate =  0.0;  CellWidth = 0.0015;}
  under ControlPointX[5] {Coordinate =  0.05; CellWidth = 0.005;}

  under ControlPointY[0] {Coordinate =  0.0;  CellWidth  = 0.0015;}
  under ControlPointY[1] {Coordinate =  0.25; CellWidth  = 0.0015;}
  under ControlPointY[2] {Coordinate =  1.0;  CellWidth  = 0.005;}
  under ControlPointY[3] {Coordinate =  2.0;  CellWidth  = 0.06;}

  NumberOfCellsZ = 1;
}


under EmbeddedBoundaryMethod {
  RiemannNormal = Mesh; //Average;
  ReconstructionAtInterface = Constant;
  under EmbeddedSurfaces {
    under Surface[0] {
      MeshFile = "wall.top";
      SurfaceThickness = 1.0e-8;
      GaussQuadrature = OnePoint;
      GaussPointsLofting = 0.0;
      InternalPressure = 1.0e5;
      //ConstantReconstructionDepth = 0.003;
    }
  }
}


under Equations {
  under Material[0] { //water, Johnsen and Colonius (2008)
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 6.59;
      PressureConstant = 4.10e8;
    }
  }

  under Material[1] { //air
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
    }
  }
}

under InitialCondition {

  UserDefinedStateCalculator = "IC/UserDefinedState.so";

  under GeometricEntities {
    under Sphere[0] {
      Center_x = -0.1;
      Center_y = 0.0;
      Center_z = 0.0;
      Radius = 0.05;
      under InitialState {
        MaterialID = 1;
        Density = 1.225e-6;
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.01e5;
      }
    }
  }
}

under BoundaryConditions {
  under Farfield {
    Density = 0.0010; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.01e5;
  }
}

under Space {
  under NavierStokes {
    Flux = LocalLaxFriedrichs;
    under Reconstruction {
      Type = Linear;    
      VariableType = ConservativeCharacteristic;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.1;
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    Solver = FiniteDifference;

    //Bandwidth = 6;

    BoundaryConditionX0   = LinearExtrapolation;
    BoundaryConditionXmax = LinearExtrapolation;
    BoundaryConditionY0   = ZeroNeumann;
    BoundaryConditionYmax = LinearExtrapolation;
    BoundaryConditionZ0   = ZeroNeumann;
    BoundaryConditionZmax = ZeroNeumann;

    under Reinitialization {
      Frequency = 4000;
      MaxIts = 800;
      FirstLayerTreatment = Fixed; //IterativelyConstrainedMethod2;
      CFL = 0.8;
      ConvergenceTolerance = 1.0e-5;
    }
  }
}

under MultiPhase {
  Flux = Numerical;
  ReconstructionAtInterface = Constant;
  PhaseChange = RiemannSolution;
}

under Time {
  Type = Explicit;
  MaxTime = 6.0e-7;
  CFL = 0.4;
  under Explicit {
    Type = RungeKutta2; 
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 1.0e-8;
  Velocity = On;
  LevelSet0 = On;
  MaterialID = On;
  Pressure = On;
  Density = On;

  under Probes {
    TimeInterval = 2.5e-10;
    //Frequency = 100;
    Pressure = "pressure_probes.txt";
    under Node[0] {X = -0.005; Y = 0.0;  Z = 0.0;}
    under Node[1] {X = -0.005; Y = 0.05; Z = 0.0;}
    under Node[2] {X = -0.005; Y = 0.1;  Z = 0.0;}
  }

  under LinePlot[0] {
    FileName = "line";
    NumberOfPoints = 1000;
    TimeInterval = 1.0e-7;
    X0   = -4.0;  Y0   = 0.0;  Z0   = 0.0;
    Xmax = 0.05;  Ymax = 0.0;  Zmax = 0.0;
  }

  MeshInformation = "meshinfo.txt";

  VerboseScreenOutput = Medium;
}

under TerminalVisualization {
  Plane = XY;
  Variable = MaterialID;
  ClockTimeInterval = 60; //seconds (computation time)
  HorizontalMin = -0.3;
  HorizontalMax = 0.05;
  VerticalMin = 0.0;
  VerticalMax = 0.2;
}

