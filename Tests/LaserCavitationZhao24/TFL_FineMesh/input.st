//units: g, mm, s, K

under Mesh {

  Type = Cylindrical;
  X0   = -6.0;
  Xmax = 6.0;
  Y0   = 0.0;
  Ymax = 6.0;
  Z0   = -0.00125;
  Zmax = 0.00125;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate = -6.0; CellWidth  = 0.2;}
  under ControlPointX[1] {Coordinate = -4.0; CellWidth  = 0.10;}
  under ControlPointX[2] {Coordinate = -2.0; CellWidth  = 0.025;}
  under ControlPointX[3] {Coordinate = -1.0; CellWidth  = 0.0075;}
  under ControlPointX[4] {Coordinate = -0.5; CellWidth  = 0.0025;}
  under ControlPointX[5] {Coordinate =  0.3; CellWidth  = 0.0025;}
  under ControlPointX[6] {Coordinate =  1.0; CellWidth  = 0.0075;}
  under ControlPointX[7] {Coordinate =  2.0; CellWidth  = 0.025;}
  under ControlPointX[8] {Coordinate =  4.0; CellWidth  = 0.10;}
  under ControlPointX[9] {Coordinate =  6.0; CellWidth  = 0.2;}

  under ControlPointY[0] {Coordinate =  0.0; CellWidth  = 0.0025;}
  under ControlPointY[1] {Coordinate =  0.4; CellWidth  = 0.0025;}
  under ControlPointY[2] {Coordinate =  1.0; CellWidth  = 0.0075;}
  under ControlPointY[3] {Coordinate =  2.0; CellWidth  = 0.025;}
  under ControlPointY[4] {Coordinate =  4.0; CellWidth  = 0.10;}
  under ControlPointY[5] {Coordinate =  6.0; CellWidth  = 0.2;}

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
    EquationOfState = NobleAbelStiffenedGas;
    under NobleAbelStiffenedGasModel {
      SpecificHeatRatio = 2.057;
      PressureConstant = 1.066e+9;
      VolumeConstant = 0.0;
      EnergyConstant = -1994.674e+9; //mm^2/s^2

      SpecificHeatAtConstantVolume = 3.449e+9;  //mm2/(s2.K)

      DensityCutOff  = 1.0e-12;
      PressureCutOff = -3.0e8;
    }

    under HeatDiffusionModel {
      Type = Constant;
      Conductivity = 0.5576e+6;
    }
  }

  under Material[1] { //water vapor
    EquationOfState = NobleAbelStiffenedGas;
    under NobleAbelStiffenedGasModel {
      SpecificHeatRatio = 1.327;
      PressureConstant = 0.0;
      VolumeConstant = 0.0;
      EnergyConstant = 1995e+9;

      SpecificHeatAtConstantVolume = 1.2e+9;

      DensityCutOff  = 1.0e-13;
      PressureCutOff = 1.0e-13;  //KW(12/2025): Try to avoid this cutoff or make it negative
    }

    under HeatDiffusionModel {
      Type = Constant;
      Conductivity = 24.57e+3;
    }
  }

  under MaterialTransition[0] { //vaporization
    FromMaterialID = 0;
    ToMaterialID = 1;
    TemperatureUpperbound = 373.15;  //assumed vaporization temperature
    LatentHeat = 2.2564e12;  //mm2/s2
  }
}

under InitialCondition {

}

under BoundaryConditions {
  under Farfield {
    Density = 0.0010; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}

under Laser {

  //SourcePower = 2.6578e12; //g.mm^2/(s^3), equiv. nanowatt 

  SourcePowerTimeHistory = "laser_power.txt";

  SourceDistribution = Gaussian;
  SourceCenterX = -0.5;
  SourceCenterY = 0.0;
  SourceCenterZ = 0.0;
  SourceBeamWaist = 0.05;
  SourceRadius = 0.1825;
  DirectionX = 1.0;
  DirectionY = 0.0;
  DirectionZ = 0.0;
  FocusingAngle = -9.78;
  Range = 5.0;

  SourceDepth = 0.025;

  //Parallelization = Original;

  Alpha = 1.0;
  RelaxationCoefficient = 1.0;
  MaxIts = 500;
  ConvergenceTolerance = 1.0e-6;

  under AbsorptionCoefficient[0] {
    MaterialID = 0;
    ReferenceCoefficient = 14; // 1/mm
  }
  under AbsorptionCoefficient[1] {
    MaterialID = 1;
    ReferenceCoefficient = 0.0001;
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
      Frequency = 1;
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

  LatentHeatTransfer = On;
}

under Time {
  Type = Explicit;
  MaxTime = 8.0e-5;
  CFL = 0.25;
  under Explicit {
    Type = RungeKutta2; 
  }
}

under ExactRiemannSolution {
  MaxIts = 200;
  MaxItsShock = 200;
  IntegrationSteps = 1000;
  Tolerance = 1.0e-9;
  ToleranceShock = 1.0e-9;
  ToleranceRarefaction = 1.0e-10;
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
*/
  MeshInformation = "meshinfo.txt";

  VerboseScreenOutput = Low;
}
