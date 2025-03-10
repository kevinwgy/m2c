//Merging and separation of two circles
//Parameters based on Section 4.2.3 of Hartmann, Meinke, Schroder 2008 (Hardcoded in M2C w/ macro LEVELSET_TEST=3)

under Mesh {
  X0   = -5.0;
  Xmax = 5.0;
  Y0   = -5.0;
  Ymax = 5.0;
  Z0   = -0.04;
  Zmax = 0.04;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Farfield;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  NumberOfCellsX = 128;
  NumberOfCellsY = 128;
  NumberOfCellsZ = 1;
}

//NOT REALLY USED
under Equations {
  under Material[0] { //air
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
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

}

under BoundaryConditions {
  under Inlet {
    Density = 0.0010; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}

under Space {
  under NavierStokes {
    Flux = HLLC;
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

    under Reinitialization {
      Frequency = 60;
      MaxIts = 500;
      ConvergenceTolerance = 1.0e-5;
      CFL = 0.8;
    }

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
  PhaseChange = RiemannSolution;
}

under Time {
  Type = Explicit;
  MaxTime = 3.2;
  TimeStep = 0.005; 
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 0.2;
  //Velocity = On;
  LevelSet0 = On;
  //MaterialID = On;

  //VerboseScreenOutput = On;
}
