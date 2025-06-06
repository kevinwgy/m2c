//units: g, mm, s, K

under Mesh {

  Type = Cylindrical;
  X0   = -2.0;
  Xmax = 0.9;
  Y0   = 0.0;
  Ymax = 1.0;
  Z0   = -0.001;
  Zmax = 0.001;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = SlipWall;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = SlipWall;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate = -2.0; CellWidth  = 0.02;}
  under ControlPointX[1] {Coordinate = -1.0; CellWidth  = 0.005;}
  under ControlPointX[2] {Coordinate = -0.5; CellWidth  = 0.0025;}
  under ControlPointX[3] {Coordinate =  0.7; CellWidth  = 0.0025;}
  under ControlPointX[4] {Coordinate =  0.9; CellWidth  = 0.004;}

  under ControlPointY[0] {Coordinate =  0.0; CellWidth  = 0.0025;}
  under ControlPointY[1] {Coordinate =  0.5; CellWidth  = 0.0025;}
  under ControlPointY[2] {Coordinate =  1.0; CellWidth  = 0.005;}

  NumberOfCellsZ = 1;
}


under EmbeddedBoundaryMethod {
  RiemannNormal = Mesh; //Average;
  ReconstructionAtInterface = Constant;
  under EmbeddedSurfaces {
    under Surface[0] {
      MeshFile = "laser_probe.top";
      SurfaceThickness = 1.0e-8;
      GaussQuadrature = OnePoint;
      GaussPointsLofting = 0.0;
      InternalPressure = 1.0e5;
      ConstantReconstructionDepth = 0.003;

      under Output {
        //Frequency = 1;
        TimeInterval = 2.0e-7;
        Prefix = "results/";
        Mesh = "surface0.top";
        Displacement = "surface0_disp.txt";
        Solution = "surface0_load.txt";
        ContactSurfaceOutput = "surface0_contact.top";
      }
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

  UserDefinedStateCalculator = "IC/libUserDefinedState.so";

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

    Bandwidth = 6;

    BoundaryConditionX0   = LinearExtrapolation;
    BoundaryConditionXmax = LinearExtrapolation;
    BoundaryConditionY0   = ZeroNeumann;
    BoundaryConditionYmax = LinearExtrapolation;
    BoundaryConditionZ0   = ZeroNeumann;
    BoundaryConditionZmax = ZeroNeumann;

    under Reinitialization {
      Frequency = 2;
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
  //MaxTime = 1.0e-6;
  MaxTime = 60.0e-6;
  CFL = 0.4;
  under Explicit {
    Type = RungeKutta2; 
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 2.0e-7;
  Velocity = On;
  LevelSet0 = On;
  MaterialID = On;
  LaserRadiance = On;
  Temperature = On;
  Pressure = On;
  Density = On;
/*
  under Probes {
    //TimeInterval = 2.0e-6;
    Frequency = 1;
    Pressure = "pressure_probes.txt";
    VelocityX = "vx_probes.txt";
    VelocityY = "vy_probes.txt";
    VelocityZ = "vz_probes.txt";
    Density  = "density_probes.txt";
    Temperature = "temperature_probes.txt";
    under Node[0] {X = 0.0; Y = 0.0; Z = 0.0;}
    under Node[1] {X = 0.5; Y = 0.005; Z = 0.0;}
  }
  
  under EnergyIntegration {
    //TimeInterval = 2.0e-6;
    Frequency = 1;
    Volume = "volume_integration.txt";
    Mass = "mass_integration.txt";
    TotalEnergy = "E_integration.txt";
    TotalEnthalpy = "H_integration.txt";
    KineticEnergy  = "kinetic_integration.txt";
    InternalEnergy = "internal_integration.txt";
    PotentialEnergy = "potential_integration.txt";

    Xmin = -12.0;
    Xmax = 12.0;
    Ymin = 0.0;
    Ymax = 12.0;
    Zmin = -1.0;
    Zmax = 1.0;
  }
*/

  MeshInformation = "meshinfo.txt";

  VerboseScreenOutput = Low;
}
