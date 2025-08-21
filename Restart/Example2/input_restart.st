//units: g, mm, s, K

under Mesh {

  Type = Cylindrical;
  X0   = -3.0;
  Xmax = 1.2;
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

  under ControlPointX[0] {Coordinate = -3.0; CellWidth  = 0.02;}
  under ControlPointX[1] {Coordinate = -2.0; CellWidth  = 0.02;}
  under ControlPointX[2] {Coordinate = -1.0; CellWidth  = 0.005;}
  under ControlPointX[3] {Coordinate = -0.5; CellWidth  = 0.0028;}
  under ControlPointX[4] {Coordinate =  0.7; CellWidth  = 0.0028;}
  under ControlPointX[5] {Coordinate =  0.9; CellWidth  = 0.004;}
  under ControlPointX[6] {Coordinate =  1.2; CellWidth  = 0.008;}

  under ControlPointY[0] {Coordinate =  0.0; CellWidth  = 0.0028;}
  under ControlPointY[1] {Coordinate =  0.5; CellWidth  = 0.0028;}
  under ControlPointY[2] {Coordinate =  1.0; CellWidth  = 0.005;}
  under ControlPointY[3] {Coordinate =  2.0; CellWidth  = 0.02;}

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
    under Surface[1] {
      MeshFile = "crater.top";
      SurfaceThickness = 1.0e-8;
      GaussQuadrature = OnePoint;
      GaussPointsLofting = 0.0;
      InternalPressure = 1.0e5;
      ConstantReconstructionDepth = 0.003;

      under Output {
        //Frequency = 1;
        TimeInterval = 2.0e-7;
        Prefix = "results/";
        Mesh = "surface1.top";
        Displacement = "surface1_disp.txt";
        Solution = "surface1_load.txt";
        ContactSurfaceOutput = "surface1_contact.top";
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
      PressureCutOff = 1.0e-13;
    }

    under HeatDiffusionModel {
      Type = Constant;
      Conductivity = 24.57e+3;
    }
  }

/*
  under MaterialTransition[0] { //vaporization
    FromMaterialID = 0;
    ToMaterialID = 1;
    TemperatureUpperbound = 373.15;  //assumed vaporization temperature
    LatentHeat = 2.2564e12;  //mm2/s2
  }
*/
}

under InitialCondition {
  UserDefinedStateCalculator = "IC/libUserDefinedState.so";

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

/*
under Laser {

  //SourcePower = 2.6578e12; //g.mm^2/(s^3), equiv. nanowatt 

  SourcePowerTimeHistory = "laser_power.txt";

  SourceDistribution = Gaussian;
  SourceCenterX = -0.5;
  SourceCenterY = 0.0;
  SourceCenterZ = 0.0;
  SourceBeamWaist = 0.165;
  SourceRadius = 0.1825;
  DirectionX = 1.0;
  DirectionY = 0.0;
  DirectionZ = 0.0;
  FocusingAngle = -7.535;
  Range = 1.1;  //cannot reach the crater wall

  SourceDepth = 0.025;


  Alpha = 1.0;
  RelaxationCoefficient = 1.0;
  MaxIts = 500;
  ConvergenceTolerance = 1.0e-6;

  under AbsorptionCoefficient[0] {
    MaterialID = 0;
    ReferenceCoefficient = 2.4202; // 1/mm
  }
  under AbsorptionCoefficient[1] {
    MaterialID = 1;
    ReferenceCoefficient = 0.0001;
  }
}
*/

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

  LatentHeatTransfer = On;
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

  MeshInformation = "meshinfo.txt";

  VerboseScreenOutput = Low;
}
