//Units: g, mm, s, Pa. Detonation problem 

under Mesh {       
  Type = Cylindrical;
  X0   =   0.0;
  Xmax =   200;
  Y0   =   0.0;
  Ymax =   200;
  Z0   =  -0.5;
  Zmax =   0.5;


  BoundaryConditionX0   = Symmetry; 
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  NumberOfCellsX = 500;
  NumberOfCellsY = 500;
  NumberOfCellsZ = 1; 

  // Uncomment the following for a non-uniform mesh.
/*
  under ControlPointX[0] {Coordinate =    0.0;   CellWidth  =  0.2;}
  under ControlPointX[1] {Coordinate =   50.0;   CellWidth  =  0.2;}
  under ControlPointX[2] {Coordinate =   100.0;  CellWidth  = 1.0;}
  under ControlPointX[3] {Coordinate =   200.0;  CellWidth  = 5.0;}

  under ControlPointY[0] {Coordinate =    0.0;   CellWidth  =  0.2;}
  under ControlPointY[1] {Coordinate =   50.0;   CellWidth  =  0.2;}
  under ControlPointY[2] {Coordinate =   100.0;  CellWidth  = 1.0;}
  under ControlPointY[3] {Coordinate =   200.0;  CellWidth  = 5.0;}
*/
}

under Equations {
  under Material[0] { // air
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0;
    }
    PressureCutOff = 1.0e-2;
    DensityCutOff  = 1.0e-11;
  }

  under Material[1] {
    EquationOfState = JonesWilkinsLee;
    under JonesWilkinsLeeModel {
      Omega = 0.28;
      A1    = 3.712e11;
      A2    = 3.21e9;
      R1    = 4.15;
      R2    = 0.95;
      Rho0  = 1.63e-3;
    }
    PressureCutOff = 1.0e-2;
    DensityCutOff = 1.0e-11;
  }
}

under InitialCondition {

  UserDefinedStateCalculator = "./IC/StateCalculator.so";

  under GeometricEntities {

    under Sphere[0] {
      Center_x = 0.0;
      Center_y = 0.0;
      Center_z = 0.0;
      Radius = 32.14;
      Side = Interior;

      under InitialState {
        MaterialID = 1;
        Density = 1.63e-3; 
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1e5;
      }
    }

  }
}

under BoundaryConditions {
  under Farfield {
    MaterialID = 0;
    Density = 1.177e-6; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}

under Space {
  under NavierStokes {
    Flux = LocalLaxFriedrichs;
    under Reconstruction {
      Type = Linear;    
      VariableType = Primitive; //sensitive, but no clear winner (depends on limiter coeff.)
      //SlopeNearInterface = NonZero;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.0;
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    Solver = FiniteDifference;

    BoundaryConditionX0   = ZeroNeumann;
    BoundaryConditionXmax = LinearExtrapolation;
    BoundaryConditionY0   = ZeroNeumann;
    BoundaryConditionYmax = LinearExtrapolation;
    BoundaryConditionZ0   = ZeroNeumann;
    BoundaryConditionZmax = ZeroNeumann;

    under Reinitialization {
      MaxIts = 1000;
      Frequency = 200;
      ConvergenceTolerance = 1e-4;
      FirstLayerTreatment = IterativelyConstrainedMethod2;
    }
  }

}

under MultiPhase {
  Flux = Numerical; //LocalLaxFriedrichs;
  ReconstructionAtInterface = Constant; //Linear; //Constant;
  PhaseChange = RiemannSolution;
  ConstantReconstructionDepth = 3.0; // this seems to give the best looking result.
}

under Time {
  Type = Explicit;
  MaxTime = 3.0e-5;
  CFL    = 0.4;
  under Explicit {
    Type = RungeKutta3;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 6.0e-7;
  Density = On;
  Velocity = On;
  Pressure = On;
  LevelSet0 = On;
  InternalEnergyPerUnitMass = On;
  MaterialID = On;

  under EnergyIntegration {
    TimeInterval = 1e-7;
    TotalEnergy = "total_energy.txt";
    KineticEnergy = "kinetic_energy.txt";
    InternalEnergy = "internal_energy.txt";
  }

  VerboseScreenOutput = Low;
}
