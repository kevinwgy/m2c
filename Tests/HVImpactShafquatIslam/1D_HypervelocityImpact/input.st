//Units: g, mm, s, Pa. 
//This is a 1D hypervelocity impact test, w/ 2 material subdomains: Tantalum and Soda-Lime Glass (SLG)

under Mesh {
  Type = ThreeDimensional;
  X0   = 0.0;
  Xmax = 1.0;
  Y0   = -0.03;
  Ymax = 0.03;
  Z0   = -0.03;
  Zmax = 0.03;
  NumberOfCellsX = 1000;
  NumberOfCellsY = 1;
  NumberOfCellsZ = 1;

  BoundaryConditionX0   = Inlet;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Symmetry;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;
}

//For material properties, refer to Shafquat Islam report (04/2021)
under Equations {

  under Material[0] { //Tantalum
    EquationOfState = MieGruneisen;
    under MieGruneisenModel {
      ReferenceDensity = 16.65e-3; //g/mm3
      BulkSpeedOfSound = 3.293e6; //mm/s
      HugoniotSlope = 1.307; //non-D
      ReferenceGamma = 1.64; //non-D
      SpecificHeatAtConstantVolume = 139e6; //139 J/(K.kg) = 139e6 mm2.s2/K (Dulong-Petit law)
      ReferenceSpecificInternalEnergy = 0.0;
      ReferenceTemperature = 270.0;
    }
    DensityCutOff = 1.0e-8;
    PressureCutOff = 1.0;
    DensityUpperLimit = 64.0e-3; //g/mm3 (M-G may give negative c^2 if 1-s*eta<0, where eta = 1-rho0/rho)
   // PressureUpperLimit = 1.0e12;
    DensityPrescribedAtFailure = 16.65e-3; //density of tantalum
  }

under Material[1] { //Soda-lime glass    //KW: check the numbers
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 3.9;
      PressureConstant = 2.6208e9; //Pa
      SpecificHeatAtConstantVolume = 1156e6; //1457.0e6; // asymptotic liquid cp: Huang and Gupta
      ReferenceTemperature = 270.0; //K
      ReferenceSpecificInternalEnergy = 1.40982e12; //e = (p + gam*pc)/((gam-1)*rho), with initial rho and p
    }
    DensityCutOff = 1.0e-6;
    PressureCutOff = -2e9;
//    PressureUpperLimit = 1.0e12;
    DensityPrescribedAtFailure = 0.0025; //density of SLG
  }

}

under InitialCondition {

  under GeometricEntities {

    under Plane[0] {
      Point_x = 0.4;
      Point_y = 0.0;
      Point_z = 0.0;
      Normal_x = 1.0;
      Normal_y = 0.0;
      Normal_z = 0.0;

      under InitialState {
        MaterialID = 1;
        Density = 2.485e-3;
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
    Density = 16.65e-3; //g/mm3
    VelocityX = 3.0e6; //mm/s
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
  under Farfield {
    MaterialID = 1;
    Density = 2.485e-3;
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
      VariableType = ConservativeCharacteristic; //C-C ~ C ~ P-C > P
      SlopeNearInterface = NonZero;
      Limiter = VanAlbada; //GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.5;
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    Flux = Roe;
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
  MaxTime = 8.0e-8;
  CFL    = 0.8;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  //Frequency = 10;
  TimeInterval = 8.0e-9;
  Density = On;
  Velocity = On;
  Pressure = On;
  LevelSet0 = On;
  MaterialID = On;
  InternalEnergyPerUnitMass = On;

  under LinePlot[0] {
    FileName = "line";
    NumberOfPoints = 200;
    TimeInterval = 8.0e-10;
    X0 = 0.0025;  Xmax = 0.9975;
    Y0 = 0.0;  Ymax = 0.0;
    Z0 = 0.0;  Zmax = 0.0;
  }

  VerboseScreenOutput = High;
}
