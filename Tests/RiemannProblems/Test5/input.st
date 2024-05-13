//Units: kg, m, s, Pa. This is the problem on Page 135 of Arthur Rallu thesis (Fig. 4.15)

under Mesh {
  Type = ThreeDimensional;
  X0   = 0.0;
  Xmax = 1.0;
  Y0   = -0.03;
  Ymax = 0.03;
  Z0   = -0.03;
  Zmax = 0.03;
  NumberOfCellsX = 200;
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
      SpecificHeatRatio = 4.4;
      PressureConstant = 6.0e8;
    }
  }
}

under InitialCondition {

  under GeometricEntities {

    under Plane[0] {
      Point_x = 0.3;
      Point_y = 0.0;
      Point_z = 0.0;
      Normal_x = 1.0;
      Normal_y = 0.0;
      Normal_z = 0.0;

      under InitialState {
        MaterialID = 1;
        Density = 1000.0; 
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.0e9;
      }
    }
  }

}

under BoundaryConditions {
  under Inlet {
    MaterialID = 0;
    Density = 50.0; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
  under Inlet2 {
    MaterialID = 1;
    Density = 1000.0; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e9;
  }
}

under Space {
  under NavierStokes {
    Flux = HLLC;
    under Reconstruction {
      Type = Linear;    
      VariableType = Conservative; //ConservativeCharacteristic; //C-C = C = P-C > P
      SlopeNearInterface = NonZero;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.0;
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    under Reconstruction {
      Type = Linear;
      Limiter = None; //Beter than VanAlbada (or any other)
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
  //MaxIts = 500;
  MaxTime = 2.4e-4;
  CFL    = 0.8;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  //Frequency = 10;
  TimeInterval = 2.4e-5;
  Density = On;
  Velocity = On;
  Pressure = On;
  LevelSet0 = On;
  MaterialID = On;

  under LinePlot[0] {
    FileName = "line";
    NumberOfPoints = 200;
//    TimeInterval = 0.001;
    X0 = 0.0;  Xmax = 1.0;
    Y0 = 0.0;  Ymax = 0.0;
    Z0 = 0.0;  Zmax = 0.0;
  }

  //VerboseScreenOutput = On;
}
