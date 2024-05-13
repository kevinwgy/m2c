//Units: kg, m, s, Pa. This is the problem on Page 142 of Arthur Rallu thesis (Fig. 4.18)

under Mesh {
  Type = ThreeDimensional;
  X0   = -0.03;
  Xmax = 0.03;
  Y0   = -0.03;
  Ymax = 0.03;
  Z0   = 0.0;
  Zmax = 1.0;
  NumberOfCellsX = 1;
  NumberOfCellsY = 1;
  NumberOfCellsZ = 200;

  BoundaryConditionX0   = Symmetry;
  BoundaryConditionXmax = Symmetry;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Symmetry;
  BoundaryConditionZ0   = Inlet;
  BoundaryConditionZmax = Inlet2;
}

under Equations {
  under Material[0] {
    EquationOfState = JonesWilkinsLee;
    under JonesWilkinsLeeModel {
      Omega = 0.28;
      A1    = 3.712e11;
      A2    = 3.23e9;
      R1    = 4.15;
      R2    = 0.95;
      Rho0  = 1630.0;
    }
  }
  under Material[1] {
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 7.15;
      PressureConstant = 3.309e8;
    }
  }
}

under InitialCondition {

  under GeometricEntities {

    under Plane[0] {
      Point_x = 0.0;
      Point_y = 0.0;
      Point_z = 0.4;
      Normal_x = 0.0;
      Normal_y = 0.0;
      Normal_z = 1.0;

      under InitialState {
        MaterialID = 1;
        Density = 1000.0; 
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.0e5;
      }
    }
  }

}

under BoundaryConditions {
  under Inlet {
    MaterialID = 0;
    Density = 1630.0; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e11;
  }
  under Inlet2 {
    MaterialID = 1;
    Density = 1000.0; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}

under Space {
  under NavierStokes {
    Flux = Roe;
    under Reconstruction {
      Type = Linear;    
      VariableType = Primitive; //sensitive, but no clear winner (depends on limiter coeff.)
      SlopeNearInterface = NonZero;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.8;
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    under Reconstruction {
      Type = Linear;
      Limiter = None;
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
  MaxTime = 3.0e-5;
  CFL    = 0.8;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  //Frequency = 10;
  TimeInterval = 3.0e-6;
  Density = On;
  Velocity = On;
  Pressure = On;
  LevelSet0 = On;
  MaterialID = On;

  under LinePlot[0] {
    FileName = "line";
    NumberOfPoints = 200;
//    TimeInterval = 0.001;
    X0 = 0.0;  Xmax = 0.0;
    Y0 = 0.0;  Ymax = 0.0;
    Z0 = 0.0025;  Zmax = 0.9975;
  }

  //VerboseScreenOutput = On;
}
