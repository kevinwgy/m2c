under Mesh {
  Type = ThreeDimensional;
  X0   = -1.0;
  Xmax = 2.0;
  Y0   = -0.05;
  Ymax = 0.05;
  Z0   = -0.05;
  Zmax = 0.05;
  NumberOfCellsX = 400;
  NumberOfCellsY = 1;
  NumberOfCellsZ = 1;

  BoundaryConditionX0   = Inlet;
  BoundaryConditionXmax = Inlet2;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Symmetry;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;
}

under Equations {
  under Material[0] {
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
    }
  }
  under Material[1] {
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
    }
  }
}

under InitialCondition {
  UserDefinedStateCalculator = "IC/libUserDefinedState.so";
}

under BoundaryConditions {
  under Inlet {
    MaterialID = 0;
    Density = 1.0; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0;
  }
  under Inlet2 {
    MaterialID = 1;
    Density = 0.125; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 0.1;
  }
}

under Space {
  under NavierStokes {
    Flux = Roe;
    under Reconstruction {
      Type = Linear;    
      VariableType = Primitive; //not very sensitive
      SlopeNearInterface = NonZero;
      Limiter = VanAlbada; //GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.5;
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    under Reconstruction {
      Type = Linear;
      Limiter = None;
    }
  }

  under PrescribedSolution {
    UserDefinedSolutionCalculator = "UDSL/UserDefinedSolution.so";
  }
}

under MultiPhase {
  Flux = Numerical;
  ReconstructionAtInterface = Linear;  //slightly better than constant
  PhaseChange = RiemannSolution; //clearly better than extrapolation
}

under Time {
  Type = Explicit;
  //MaxIts = 500;
  MaxTime = 0.1;
  //MaxTime = 0.2;
  CFL    = 0.8;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  //Frequency = 10;
  TimeInterval = 0.01;
  Density = On;
  Velocity = On;
  Pressure = On;
  LevelSet0 = On;
  MaterialID = On;

  under LinePlot[0] {
    FileName = "line";
    NumberOfPoints = 200;
    Frequency = 1; //I just want to see the solution after 1 time-step...
//    TimeInterval = 0.001;
    X0 = 0.0025;  Xmax = 0.9975;
    Y0 = 0.0;  Ymax = 0.0;
    Z0 = 0.0;  Zmax = 0.0;
  }

  //VerboseScreenOutput = On;
}
