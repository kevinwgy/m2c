//Units: mm, g, s

under Mesh {
  Type = Cylindrical;
  X0   = -2000.0;
  Xmax = 2000.0;
  Y0   = 0.0;
  Ymax = 2000.0;
  Z0   = -2.5;
  Zmax = 2.5;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield2;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = SlipWall;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate =-2000.0; CellWidth  = 60.0;}
  under ControlPointX[1] {Coordinate =-1000.0; CellWidth  = 15.0;}
  under ControlPointX[2] {Coordinate = -500.0; CellWidth  = 2.0;}
  under ControlPointX[3] {Coordinate =  700.0; CellWidth  = 2.0;}
  under ControlPointX[4] {Coordinate = 1000.0; CellWidth  = 15.0;}
  under ControlPointX[5] {Coordinate = 2000.0; CellWidth  = 60.0;}

  under ControlPointY[0] {Coordinate =    0.0; CellWidth  = 2.0;}
  under ControlPointY[1] {Coordinate =  700.0; CellWidth  = 2.0;}
  under ControlPointY[2] {Coordinate = 1000.0; CellWidth  = 15.0;}
  under ControlPointY[3] {Coordinate = 2000.0; CellWidth  = 60.0;}

  NumberOfCellsZ = 1;
}

under Equations {
  under Material[0] { //Air
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0;

      // Apply a simple temperature law, not "thermodynamically consistent"
      SpecificHeatAtConstantVolume = 717.1e6;  // 717.1 J/(K.kg) = 717.1e6 mm2.s2/K
      ReferenceTemperature = 300.0; //K
      ReferenceSpecificInternalEnergy = 2.083333e11; //1.0e5/((1.4-1.0)*1.2e-6) ensures T = T0 in ref. cond.
    }
    //DensityCutOff = 1.0e-12;
    //PressureCutOff = 1.0;
    //DensityPrescribedAtFailure = 1.2e-6; //density of air
  }

  under Material[1] { //Water 

    under NobleAbelStiffenedGasModel { //Zein et al., 2013, used in Zhao et al., 2024
      SpecificHeatRatio = 2.057;
      PressureConstant = 1.066e+9;
      VolumeConstant = 0.0;
      EnergyConstant = -1994.674e+9; //mm^2/s^2

      SpecificHeatAtConstantVolume = 3.449e+9;  //mm2/(s2.K)

      //DensityCutOff  = 1.0e-12;
      //PressureCutOff = -3.0e8;
    }
  }
}


under InitialCondition {

  under GeometricEntities {
    under Plane[0] {
      Point_x = 0.0;
      Point_y = 0.0;
      Point_z = 0.0;
      Normal_x = 1.0;
      Normal_y = 0.0;
      Normal_z = 0.0;
      under InitialState {
        MaterialID = 1; //water
        Density = 0.001;
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
    Density = 1.2e-6; //density of air
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
  under Farfield2 {
    Density = 0.001; //density of water
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}

under EmbeddedBoundaryMethod {
  RiemannNormal = Mesh; //Average;
  ReconstructionAtInterface = Constant;
  under EmbeddedSurfaces {
    under Surface[0] {
      MeshFile = "impactor.top";
      SurfaceThickness = 1.0e-3;
      GaussQuadrature = OnePoint;
      GaussPointsLofting = 0.0;
      InternalPressure = 1.0e5;
      //ConstantReconstructionDepth = 20.0;

      UserDefinedDynamicsCalculator = "impactor_dynamics/UserDefinedDynamics.so";

      under Output {
        TimeInterval = 2.0e-5;
        Prefix = "results/";
        Mesh = "surface.top";
        Displacement = "surface_disp.txt";
        Solution = "surface_load.txt";
        ContactSurfaceOutput = "surface_contact.top";
      }

    }
  }
}

under Space {
  under NavierStokes {
    Flux = LocalLaxFriedrichs; //HLLC;
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

    BoundaryConditionX0   = LinearExtrapolation;
    BoundaryConditionXmax = LinearExtrapolation;
    BoundaryConditionY0   = ZeroNeumann;
    BoundaryConditionYmax = LinearExtrapolation;
    BoundaryConditionZ0   = ZeroNeumann;
    BoundaryConditionZmax = ZeroNeumann;

    //Bandwidth = 7;
    under Reinitialization {
      Frequency = 5000;
      MaxIts = 200;
    }
  }

}

under MultiPhase {
  Flux = Numerical;
  ReconstructionAtInterface = Constant;
  PhaseChange = RiemannSolution;
  RiemannNormal = Average;

  LevelSetCorrectionFrequency = 200;
  //ConstantReconstructionDepth = 10.0; //mm
}

under Time {
  Type = Explicit;
  //MaxIts = 2560;
  MaxTime = 1.2e-3;
  CFL    = 0.2;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 2.0e-5;
  Density = On;
  Velocity = On;
  Pressure = On;
  LevelSet0 = On;
  Temperature = On;
  MaterialID = On;

  VerboseScreenOutput = Medium; //Low; //High;
}

under TerminalVisualization {
  Plane = XY;
  Variable = MaterialID; //Pressure; //Velocity;
  ClockTimeInterval = 20;
  HorizontalMin = -800.0;
  HorizontalMax = 800.0;
  VerticalMin = 0.0;
  VerticalMax = 1000.0;
}

