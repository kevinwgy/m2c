//Units: g, mm, s, Pa. This is a hypervelocity impact test. The initial condition features a big jump in velocity

under Mesh {
  Type = ThreeDimensional;
  X0   = -1.0;
  Xmax = 2.0;
  Y0   = -0.03;
  Ymax = 0.03;
  Z0   = -0.03;
  Zmax = 0.03;
  NumberOfCellsX = 5000;
  NumberOfCellsY = 1;
  NumberOfCellsZ = 1;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield2;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Symmetry;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;
}

under Equations {
  

  under Material[1] { //Tantalum
    EquationOfState = ExtendedMieGruneisen;
    under ExtendedMieGruneisenModel {
      ReferenceDensity = 16.65e-3; //g/mm3
      BulkSpeedOfSound = 3.293e6; //mm/s
      HugoniotSlope = 1.307; //non-D
      ReferenceGamma = 1.64; //non-D
      SpecificHeatAtConstantVolume = 139e6; //139 J/(K.kg) = 139e6 mm2.s2/K (Dulong-Petit law)
      ReferenceSpecificInternalEnergy = 0.0;
      ReferenceTemperature = 270.0;

      //TemperatureLaw = SimplifiedCv;

      VolumetricStrainBreak = -4.0;
    }
    DensityCutOff = 1.0e-8;
    //PressureCutOff = 1.0;
    DensityUpperLimit = 64.0e-3; //g/mm3 (M-G may give negative c^2 if 1-s*eta<0, where eta = 1-rho0/rho)
    //PressureUpperLimit = 1.0e12;
    DensityPrescribedAtFailure = 16.65e-3; //density of tantalum
  }


  under Material[2] { // SLG
    EquationOfState = NobleAbelStiffenedGas;
    under NobleAbelStiffenedGasModel {
      SpecificHeatRatio = 1.2229; // gamma
      PressureConstant = 1.4136e10; // pc
      VolumeConstant = 1.46e2; // b
      EnergyConstant = 0.0; //mm^2/s^2 // q
      SpecificHeatAtConstantVolume = 2.4204e9;  //mm2/(s2.K)
      TemperatureLawConstant = 2.0613e2; //mm3/g
      DensityCutOff  = 1.0e-12;
//      PressureCutOff = -3.0e8;
    }
  }


  under Material[0] { //Ar
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.6666667;
      PressureConstant = 0;
      SpecificHeatAtConstantVolume = 312.2e6;  // 312.2 J/(K.kg) = 312.2e6 mm2.s2/K
      // e0 = 0, T0 = 0 ==> when p = 1e5, rho = 1.78e-6, T = 270 K 
    }
    DensityCutOff = 1.0e-10;
    PressureCutOff = 10.0;
//    PressureUpperLimit = 1.0e12;
    DensityPrescribedAtFailure = 1.78e-6; //density of Ar
  }
}

under InitialCondition {

  under GeometricEntities {

    under CylinderAndCone[0] { //Tantalum projectile
      Axis_x = 1.0;
      Axis_y = 0.0;
      Axis_z = 0.0;
      BaseCenter_x = 0.30;
      BaseCenter_y = 0.0;
      BaseCenter_z = 0.0;
      CylinderRadius = 30.0;
      CylinderHeight =0.20;
      under InitialState {
        MaterialID = 1;
        Density = 16.65e-3;
        VelocityX = 5.0e6;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.0e5;
      }
    }

    under CylinderAndCone[1] { //SLG target
      Axis_x = 1.0;
      Axis_y = 0.0;
      Axis_z = 0.0;
      BaseCenter_x = 0.50;
      BaseCenter_y = 0.0;
      BaseCenter_z = 0.0;
      CylinderRadius = 30.0;
      CylinderHeight =0.25;
      under InitialState {
        MaterialID = 2;
        Density = 2.485e-3;
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.0e5;
      }
    }

    under CylinderAndCone[2] { //Outlet Ar
      Axis_x = 1.0;
      Axis_y = 0.0;
      Axis_z = 0.0;
      BaseCenter_x = 0.75;
      BaseCenter_y = 0.0;
      BaseCenter_z = 0.0;
      CylinderRadius = 30.0;
      CylinderHeight =2.0;
      under InitialState {
        MaterialID = 0;
        Density = 1.78e-6;
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.0e5;
      }
    }

  }

}

under BoundaryConditions {
  under Farfield {
    MaterialID = 0;
    Density = 1.78e-6; //g/mm3
    VelocityX = 5.0e6; //mm/s
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
  under Farfield2 {
    MaterialID = 0;
    Density = 1.78e-6;
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}

under Space {
  under NavierStokes {
    Flux = LocalLaxFriedrichs; //HLLC; //Godunov; //LocalLaxFriedrichs;
    under Reconstruction {
      Type = Linear; //Constant;
      VariableType = ConservativeCharacteristic; //C-C ~ C ~ P-C > P
      SlopeNearInterface = NonZero;
      Limiter = GeneralizedMinMod; //VanAlbada; //GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.1;
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    under Reconstruction {
      Type = Linear;
      Limiter = None;
    }
  }
  under LevelSet[1] {
    MaterialID = 2;
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
  CFL    = 0.1;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  //Frequency = 10;
  TimeInterval = 1.0e-8;
  Density = On;
  Velocity = On;
  Pressure = On;
  LevelSet0 = On;
  MaterialID = On;

  under LinePlot[0] {
    FileName = "line";
    NumberOfPoints = 10000;
    TimeInterval = 1e-10;
    X0 = -0.9975;  Xmax = 1.9975;
    Y0 = 0.0;  Ymax = 0.0;
    Z0 = 0.0;  Zmax = 0.0;
  }

  //VerboseScreenOutput = On;
}
