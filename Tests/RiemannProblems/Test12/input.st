// This is a test case featuring a solid flying away from air at very high velocity (Actually , if the velocity is above 1.8e6 mm/s, the problem has no solution.
// units: g, mm, s

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
  under Material[0] { //air
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
      DensityCutOff = 1.0e-13;
      PressureCutOff = 1.0e-3;
    }
  }
  under Material[1] { //SLG (Starphire) 
    EquationOfState = MieGruneisen;
    under MieGruneisenModel {
      ReferenceDensity = 2.204e-3; //g/mm3   CORRECT???
      BulkSpeedOfSound = 2.22e6; //mm/s
      HugoniotSlope = 1.61; //non-D
      ReferenceGamma = 0.65; //non-D   CORRECT???
      PressureCutOff = 1.0e-3;
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
        Density = 2.204e-3; //g/mm3
        VelocityX = 1.5e6;
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
    Density = 1.2e-6; //g/mm3
    VelocityX = 0.0; //mm/s
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
  under Inlet2 {
    MaterialID = 1;
    Density = 2.204e-3; //g/mm3
    VelocityX = 1.5e6;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}

under Space {
  under NavierStokes {
    Flux = LocalLaxFriedrichs; //only LLF works
    under Reconstruction {
      Type = Linear;    
      VariableType = PrimitiveCharacteristic; //PrimitiveCharacteristic ~ C-C, slighter better than P. C does not work
      SlopeNearInterface = NonZero;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.3;
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
  MaxTime = 1.5e-7;
  CFL    = 0.5;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  //Frequency = 10;
  TimeInterval = 2.5e-8;
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
