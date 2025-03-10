//Vortex deformation of a circle 
//Parameters based on Section 5.2 of Hartmann, Meinke, Schroder 2010 (Hardcoded in M2C w/ macro LEVELSET_TEST=2)

under Mesh {
  X0   = 0.0;
  Xmax = 1.0;
  Y0   = 0.0;
  Ymax = 1.0;
  Z0   = -2e-3;
  Zmax = 2e-3;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Farfield;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  NumberOfCellsX = 256;
  NumberOfCellsY = 256;
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

    //Bandwidth = 6;

    under Reinitialization {
      Frequency = 10;
      MaxIts = 500;
      CFL = 0.8;
      ConvergenceTolerance = 1.0e-5;
      FirstLayerTreatment = IterativelyConstrainedMethod2;
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
  MaxTime = 8.0;
  TimeStep = 3.0e-3; 
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 0.4;
  //Velocity = On;
  LevelSet0 = On;
  //MaterialID = On;

  //VerboseScreenOutput = On;
}
