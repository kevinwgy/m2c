//Units: non-dimensional. This is the problem on Page 133 of Arthur Rallu thesis (Fig. 4.14)

under Mesh {
  Type = ThreeDimensional;
  X0   = 0.0;
  Xmax = 1.0;
  Y0   = -0.02;
  Ymax = 0.02;
  Z0   = -0.02;
  Zmax = 0.02;
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
      SpecificHeatRatio = 1.6666667;
      PressureConstant = 0.0;
    }
  }
  under Material[1] {
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.2;
      PressureConstant = 0.0;
    }
  }
}

under InitialCondition {

  under GeometricEntities {

    under Plane[0] {
      Point_x = 0.2;
      Point_y = 0.0;
      Point_z = 0.0;
      Normal_x = 1.0;
      Normal_y = 0.0;
      Normal_z = 0.0;

      under InitialState {
        MaterialID = 1;
        Density = 1.0; 
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.0;
      }
    }
  }

}

under BoundaryConditions {
  under Inlet {
    MaterialID = 0;
    Density = 3.2; 
    VelocityX = 9.43499279;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 100.0;
  }
  under Inlet2 {
    MaterialID = 1;
    Density = 1.0; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0;
  }
}

under Space {
  under NavierStokes {
    Flux = HLLC; //clearly better than Roe
    under Reconstruction {
      Type = Linear;    
      VariableType = Primitive; //clearly better than others
      SlopeNearInterface = NonZero;
      Limiter = VanAlbada; //GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.2;
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    Solver = FiniteVolume; //for this case, FV works better than FD (not sure why...)
    Flux = Roe;
    under Reconstruction {
      Type = Linear;
      Limiter = None; //VanAlbada;
    }
  }

}

under MultiPhase {
  Flux = Numerical;
  ReconstructionAtInterface = Constant;
  PhaseChange = Extrapolation; //clearly better than RiemannSolution
}

under Time {
  Type = Explicit;
  //MaxIts = 500;
  MaxTime = 0.06;
  CFL    = 0.8;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  //Frequency = 10;
  TimeInterval = 0.001;
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
