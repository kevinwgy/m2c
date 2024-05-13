//Units: g, mm, s, Pa. This is another hypervelocity impact test (copper-copper). The initial condition features only a big jump in velocity, while the density and pressure fields are constant

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

//For material properties, refer to Shafquat Islam report (04/2021)
under Equations {
  under Material[0] { //OFHC Copper 
    EquationOfState = MieGruneisen;
    under MieGruneisenModel {
      ReferenceDensity = 8.96e-3; //g/mm3
      BulkSpeedOfSound = 3.97e6; //mm/s
      HugoniotSlope = 1.479; //non-D
      ReferenceGamma = 2.12; //non-D
    }
  }
  under Material[1] { //OFHC Copper 
    EquationOfState = MieGruneisen;
    under MieGruneisenModel {
      ReferenceDensity = 8.96e-3; //g/mm3
      BulkSpeedOfSound = 3.97e6; //mm/s
      HugoniotSlope = 1.479; //non-D
      ReferenceGamma = 2.12; //non-D
    }
  }
}

under InitialCondition {

  under GeometricEntities {

    under Plane[0] {
      Point_x = 0.5;
      Point_y = 0.0;
      Point_z = 0.0;
      Normal_x = 1.0;
      Normal_y = 0.0;
      Normal_z = 0.0;

      under InitialState {
        MaterialID = 1;
        Density = 8.96e-3; //g/mm3
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
    Density = 8.96e-3; //g/mm3
    VelocityX = 5.0e6; //mm/s
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
  under Inlet2 {
    MaterialID = 1;
    Density = 8.96e-3; //g/mm3
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
      VariableType = ConservativeCharacteristic;
      SlopeNearInterface = NonZero;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 2.0;
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
  MaxTime = 2.5e-8;
  CFL    = 0.1;
  under Explicit {
    Type = RungeKutta3;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  //Frequency = 10;
  TimeInterval = 2.5e-9;
  Density = On;
  Velocity = On;
  Pressure = On;
  LevelSet0 = On;
  MaterialID = On;

  under LinePlot[0] {
    FileName = "line";
    NumberOfPoints = 200;
    X0 = 0.0025;  Xmax = 0.9975;
    Y0 = 0.0;  Ymax = 0.0;
    Z0 = 0.0;  Zmax = 0.0;
  }

  //VerboseScreenOutput = On;
}
