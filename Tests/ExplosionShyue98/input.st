// Ref: Wentao Ma et al., JCP 2023
// This is the test case in Shyue, 1998

under Mesh {
  Type = ThreeDimensional;
  X0   = 0.0;
  Xmax = 1.0;
  Y0   = 0.0;
  Ymax = 1.0;
  Z0   = -0.1;
  Zmax = 0.1;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Farfield;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  NumberOfCellsX = 500;
  NumberOfCellsY = 500;
  NumberOfCellsZ = 1;
}

under Equations {
  under Material[0] { // liquid
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 5.5;
      PressureConstant = 1.505;
    }
  }
  under Material[1] { // gas bubble 
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
    } 
  }
}


under InitialCondition {
  under GeometricEntities {
    under Sphere[0] {
      Center_x = 0.5;
      Center_y = 0.5;
      Center_z = 0.0;
      Radius = 0.2;
      under InitialState {
        MaterialID = 1;
        Density = 1.241; 
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 2.753;
      }
    }
  }
}

under BoundaryConditions {
  under Farfield {
    Density = 0.991; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 3.059e-4;
  }
}

under Space {
  under NavierStokes {
    Flux = Roe; // LocalLaxFriedrichs;
    under Reconstruction {
      Type = Linear;    
      VariableType = ConservativeCharacteristic;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.2;
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    Solver = FiniteDifference;

    BoundaryConditionX0   = LinearExtrapolation;
    BoundaryConditionXmax = LinearExtrapolation;
    BoundaryConditionY0   = LinearExtrapolation;
    BoundaryConditionYmax = LinearExtrapolation;
    BoundaryConditionZ0   = ZeroNeumann;
    BoundaryConditionZmax = ZeroNeumann;
  }
}

under MultiPhase {
  Flux = Numerical;
  ReconstructionAtInterface = Constant;
  PhaseChange = Extrapolation; //RiemannSolution;
}

under Time {
  Type = Explicit;
  MaxTime = 0.058;
  CFL = 0.4;
  under Explicit {
    Type = RungeKutta3; //ForwardEuler; RungeKutta2
  }
}

under ExactRiemannSolution {
  MaxIts = 200;
  MaxItsShock = 200;
  IntegrationSteps = 200;
  Tolerance = 1.0e-9;
  ToleranceShock = 1.0e-9;
  ToleranceRarefaction = 1.0e-9;
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 0.01;
  Velocity = On;
  LevelSet0 = On;
  MaterialID = On;
  Pressure = On;
  Density = On;

  under LinePlot[0] {
    FileName = "line";
    NumberOfPoints = 500;
    X0 = 0.0;  Xmax = 1.0;
    Y0 = 0.5;  Ymax = 0.5;
    Z0 = 0.0;  Zmax = 0.0;
  }
}

under TerminalVisualization {
  Plane = XY;
  Variable = Pressure; //Velocity;
  ClockTimeInterval = 10;
}


