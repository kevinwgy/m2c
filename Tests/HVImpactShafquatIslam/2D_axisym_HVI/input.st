//Units: mm, g, s, K, A
//This is a hypervelocity impact test, w/ 3 material subdomains: ambient fluid > Ar (PG EOS) |       projectile > Ta (MG EOS)        |       target > SLG (NASG EOS)
//The initial condition features a big jump in velocity: 6km/s projectile and  stationary target and fluid behind target

under Mesh { 
  Type = Cylindrical;
  X0   = -100.0;
  Xmax = 120.0;
  Y0   = 0.0;
  Ymax = 120.0;
  Z0   = -5.0e-2;
  Zmax = 5.0e-2;

  BoundaryConditionX0   = Inlet;
  BoundaryConditionXmax = Inlet;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Inlet;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate =-100.0; CellWidth  = 3.00;}
  under ControlPointX[1] {Coordinate =-50.0; CellWidth  = 0.75;}
  under ControlPointX[2] {Coordinate =-35.0; CellWidth  = 0.25;}
  under ControlPointX[3] {Coordinate =-20.0; CellWidth  = 0.1;}
  under ControlPointX[4] {Coordinate =  0.0; CellWidth  = 0.1;}
  under ControlPointX[5] {Coordinate = 20.0; CellWidth  = 0.1;}
  under ControlPointX[6] {Coordinate = 75.0; CellWidth  = 0.25;}
//  under ControlPointX[7] {Coordinate = 75.0; CellWidth  = 0.75;}
  under ControlPointX[8] {Coordinate = 120.0; CellWidth  = 3.00;}

  under ControlPointY[0] {Coordinate =  0.0; CellWidth  = 0.1;}
  under ControlPointY[1] {Coordinate = 20.0; CellWidth  = 0.1;}
  under ControlPointY[2] {Coordinate = 75.0; CellWidth  = 0.25;}
//  under ControlPointY[3] {Coordinate = 75.0; CellWidth  = 0.75;}
  under ControlPointY[4] {Coordinate = 120.0; CellWidth  = 3.00;}

  NumberOfCellsZ = 1;
}

under Equations {
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

  under Material[1] { // SLG
    EquationOfState = NobleAbelStiffenedGas;
    under NobleAbelStiffenedGasModel {
      SpecificHeatRatio = 1.2229; // gamma
      PressureConstant = 1.4136e10; // pc
      VolumeConstant = 1.46e2; // b
      EnergyConstant = 0.0; //mm^2/s^2 // q
      SpecificHeatAtConstantVolume = 2.4204e9;  //mm2/(s2.K)
      TemperatureLawConstant = 2.0613e2; //mm3/g
      DensityCutOff  = 1.0e-12;
      PressureCutOff = -3.0e8;
    }
  }

  under Material[2] { //Tantalum
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
    PressureCutOff = 1.0;
    DensityUpperLimit = 64.0e-3; //g/mm3 (M-G may give negative c^2 if 1-s*eta<0, where eta = 1-rho0/rho)
   // PressureUpperLimit = 1.0e12;
    DensityPrescribedAtFailure = 16.65e-3; //density of tantalum
  }
}

under Ionization {

  PlanckConstant = 6.62607004e-25;  // g.mm2/s
  ElectronCharge = 1.60217662e-19;  // A.s  (Coulumbs)
  ElectronMass   = 9.10938356e-28;  // g
  BoltzmannConstant = 1.38064852e-14;  // g.mm2/(s2.K)

  under Material[0] {
    Type = IdealSahaEquation;
    PartitionFunctionEvaluation = OnTheFly; //CubicSplineInterpolation; //LinearInterpolation; //OnTheFly;
    MaxIts = 200;
    ConvergenceTolerance = 1.0e-10;

    under Element[0] { //Ar
      MolarFraction = 1.0;
      MolarMass = 39.948; //g/mol
      AtomicNumber = 18;
      MaxChargeNumber = 18;
      IonizationEnergyFile = "AtomicData_g_mm_s_K_A/I_18.txt";
      ExcitationEnergyFilesPrefix = "AtomicData_g_mm_s_K_A/E_18_";
      ExcitationEnergyFilesSuffix = ".txt";
      DegeneracyFilesPrefix  = "AtomicData_g_mm_s_K_A/g_18_";
      DegeneracyFilesSuffix  = ".txt";
    }
  }

  under Material[1] {
    Type = NonIdealSahaEquation; //IdealSahaEquation;
    DepressionModel = Ebeling; //Griem; //Ebeling; //None;
    PartitionFunctionEvaluation = OnTheFly; //CubicSplineInterpolation; //LinearInterpolation;
    MaxIts = 200;
    ConvergenceTolerance = 1.0e-5;

    under Element[0] { //O
      MolarFraction = 0.601;
      MolarMass = 15.999; //g/mol
      AtomicNumber = 8;
      MaxChargeNumber = 8;
      IonizationEnergyFile = "AtomicData_g_mm_s_K_A/I_8.txt";
      ExcitationEnergyFilesPrefix = "AtomicData_g_mm_s_K_A/E_8_";
      ExcitationEnergyFilesSuffix = ".txt";
      DegeneracyFilesPrefix  = "AtomicData_g_mm_s_K_A/g_8_";
      DegeneracyFilesSuffix  = ".txt";
    }

    under Element[1] { //Si
      MolarFraction = 0.249;
      MolarMass = 20.086; //g/mol
      AtomicNumber = 14;
      MaxChargeNumber = 14;
      IonizationEnergyFile = "AtomicData_g_mm_s_K_A/I_14.txt";
      ExcitationEnergyFilesPrefix = "AtomicData_g_mm_s_K_A/E_14_";
      ExcitationEnergyFilesSuffix = ".txt";
      DegeneracyFilesPrefix  = "AtomicData_g_mm_s_K_A/g_14_";
      DegeneracyFilesSuffix  = ".txt";    
    }

    under Element[2] { //Na
      MolarFraction = 0.093;
      MolarMass = 22.99; //g/mol
      AtomicNumber = 11;
      MaxChargeNumber = 11;
      IonizationEnergyFile = "AtomicData_g_mm_s_K_A/I_11.txt";
      ExcitationEnergyFilesPrefix = "AtomicData_g_mm_s_K_A/E_11_";
      ExcitationEnergyFilesSuffix = ".txt";
      DegeneracyFilesPrefix  = "AtomicData_g_mm_s_K_A/g_11_";
      DegeneracyFilesSuffix  = ".txt";     
    }
    under Element[3] { //Ca
      MolarFraction = 0.055;
      MolarMass = 40.078; //g/mol
      AtomicNumber = 20;
      MaxChargeNumber = 20;
      IonizationEnergyFile = "AtomicData_g_mm_s_K_A/I_20.txt";
      ExcitationEnergyFilesPrefix = "AtomicData_g_mm_s_K_A/E_20_";
      ExcitationEnergyFilesSuffix = ".txt";
      DegeneracyFilesPrefix  = "AtomicData_g_mm_s_K_A/g_20_";
      DegeneracyFilesSuffix  = ".txt";     
    }
    under Element[4] { //Mg
      MolarFraction = 0.002;
      MolarMass = 24.305; //g/mol
      AtomicNumber = 12;
      MaxChargeNumber = 12;
      IonizationEnergyFile = "AtomicData_g_mm_s_K_A/I_12.txt";
      ExcitationEnergyFilesPrefix = "AtomicData_g_mm_s_K_A/E_12_";
      ExcitationEnergyFilesSuffix = ".txt";
      DegeneracyFilesPrefix  = "AtomicData_g_mm_s_K_A/g_12_";
      DegeneracyFilesSuffix  = ".txt";     
    }
  }

  under Material[2] {
    Type = NonIdealSahaEquation; //IdealSahaEquation;
    DepressionModel = Ebeling; //Ebeling; //Griem; //Ebeling; //None;
    PartitionFunctionEvaluation = OnTheFly; //CubicSplineInterpolation; //LinearInterpolation;
    MaxIts = 200;
    ConvergenceTolerance = 1.0e-5;

    under Element[0] { //Ta
      MolarFraction = 1;
      MolarMass = 180.95; //g/mol
      AtomicNumber = 73;
      MaxChargeNumber = 73;
      IonizationEnergyFile = "AtomicData_g_mm_s_K_A/I_73.txt";
      ExcitationEnergyFilesPrefix = "AtomicData_g_mm_s_K_A/E_73_";
      ExcitationEnergyFilesSuffix = ".txt";
      DegeneracyFilesPrefix  = "AtomicData_g_mm_s_K_A/g_73_";
      DegeneracyFilesSuffix  = ".txt";
    }
  }  
}

under InitialCondition {

  under GeometricEntities {

    under CylinderWithSphericalCaps[0] { //tantalum projectile
      Axis_x = 1.0;
      Axis_y = 0.0;
      Axis_z = 0.0;
      BaseCenter_x = -113.105;
      BaseCenter_y = 0.0;
      BaseCenter_z = 0.0;
      CylinderRadius = 5.025; //make it slightly bigger to "cover" the user_specified data
      CylinderHeight = 106.08;
      FrontSphericalCap = On;
      BackSphericalCap = On;
      under InitialState {
        MaterialID = 2;
        Density = 16.65e-3;
        VelocityX = 6.0e6; //mm/s
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.0e5;
      }
    }

    under CylinderAndCone[1] { //SLG target
      Axis_x = 1.0;
      Axis_y = 0.0;
      Axis_z = 0.0;
      BaseCenter_x = 0.0;
      BaseCenter_y = 0.0;
      BaseCenter_z = 0.0;
      CylinderRadius = 30.0;
      CylinderHeight = 20.0;
      under InitialState {
        MaterialID = 1;
        Density = 0.0025;
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.0e5;
      }
    }
  }

//  UserDataFile = "ic.txt";
  ApplyUserDataBeforeGeometricEntities = Yes;
}

under BoundaryConditions {
  under Inlet {
    Density = 1.78e-6; //density of Ar
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }

  under GeometricEntities2D {
    under Disk[0] {
      Center_x = -100.0;
      Center_y = 0.0;
      Center_z = 0.0;
      Normal_x = 1.0;
      Normal_y = 0.0;
      Normal_z = 0.0;
      Radius   = 5.025; //must match the radius of the projectile

      under BoundaryState { //the state of the projectile
        MaterialID = 2;
        Density = 16.65e-3;
        VelocityX = 6.0e6; //mm/s
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1.0e5;
      }
    }
  }
}

under Space {
  under NavierStokes {
    Flux = LocalLaxFriedrichs; //HLLC;
    under Reconstruction {
      Type = Linear; //Constant; //Linear;   // KW: we want linear recon. Other parameters can be adjusted. 
      VariableType = ConservativeCharacteristic;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.1;
    }
/*
    under Smoothing {
      Type = Gaussian;
      TimeInterval = 1.0e-7; 
      SigmaFactor = 0.5;
    }
*/
  }

  under LevelSet[0] {   
    MaterialID = 1;
    Solver = FiniteDifference;

    Bandwidth = 15; //7

    BoundaryConditionX0   = LinearExtrapolation;
    BoundaryConditionXmax = LinearExtrapolation;
    BoundaryConditionY0   = ZeroNeumann;
    BoundaryConditionYmax = LinearExtrapolation;
    BoundaryConditionZ0   = ZeroNeumann;
    BoundaryConditionZmax = ZeroNeumann;

    under Reinitialization {
      Frequency = 2;
      MaxIts = 200;
    }
  }

  under LevelSet[1] {   
    MaterialID = 2;
    Solver = FiniteDifference;

    Bandwidth = 15; //7

    BoundaryConditionX0   = ZeroNeumann; //NOTE: the projectile touches the left boundary!
    BoundaryConditionXmax = LinearExtrapolation;
    BoundaryConditionY0   = ZeroNeumann;
    BoundaryConditionYmax = LinearExtrapolation;
    BoundaryConditionZ0   = ZeroNeumann;
    BoundaryConditionZmax = ZeroNeumann;

    under Reinitialization {
      Frequency = 2;
      MaxIts = 200;
    }
  }

}

under MultiPhase {  //KW: feel free to change
  Flux = Numerical;
  ReconstructionAtInterface = Constant;
  PhaseChange = RiemannSolution;
  RiemannNormal = Average;

  LevelSetCorrectionFrequency = 100;

  ConstantReconstructionDepth = 0.3; //mm
}

under Time {
  Type = Explicit;
  //MaxIts = 2560;
  MaxTime = 1.2e-4;
  CFL    = 0.1;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
//  Frequency = 10;
  TimeInterval = 2.5e-8; //1e-10; //2.5e-8;
  Density = On;
  Velocity = On;
  Pressure = On;
  LevelSet0 = On;
  LevelSet1 = On;
  Temperature = On;
  DeltaInternalEnergyPerUnitMass = On;
  MaterialID = On;
  MeanCharge = On;
  ElectronDensity = On;
//  MaxChargeNumber = 2;
//  MolarFractionsElement0 = On;

  under EnergyIntegration {
    Frequency = 1;
    Volume = "Volume_integral.txt";
    Mass =  "Mass_integral.txt";
    TotalEnergy =  "TotalEnergy_integral.txt";
    TotalEnthalpy =  "TotalEnthalpy_integral.txt";
    KineticEnergy = "KineticEnergy_integral.txt";
    InternalEnergy = "InternalEnergy_integral.txt";
    PotentialEnergy = "PotentialEnergy_integral.txt";
    Xmin = -100.0;
    Xmax = 120.0;
    Ymin = 0.0;
    Ymax = 120.0;
    Zmin = -5.0e-2;
    Zmax = 5.0e-2;

  }
  under LinePlot[0] {
    FileName = "line";
    TimeInterval = 1e-8;
    NumberOfPoints = 8000;
    X0 = -5;  Xmax = 30;
    Y0 = 0.0;  Ymax = 0.0;
    Z0 = 0.0;  Zmax = 0.0;
  }


  under Probes {
    Frequency = 10;
    Pressure = "pressure_probes.txt";
    Density = "density_probes.txt";
    Temperature = "temperature_probes.txt"; 
    IonizationResult = "ionization_probes.txt";



    under Node[0] {X = -4.0; Y = 0.0; Z = 0.0;}
    under Node[1] {X = -4.0; Y = 2.0; Z = 0.0;}
    under Node[2] {X = -4.0; Y = 4.0; Z = 0.0;}
    under Node[3] {X = -4.0; Y = 6.0; Z = 0.0;}
    under Node[4] {X = -4.0; Y = 8.0; Z = 0.0;}
    
    under Node[5] {X = -2.0; Y = 0.0; Z = 0.0;}
    under Node[6] {X = -2.0; Y = 2.0; Z = 0.0;}
    under Node[7] {X = -2.0; Y = 4.0; Z = 0.0;}
    under Node[8] {X = -2.0; Y = 6.0; Z = 0.0;}
    under Node[9] {X = -2.0; Y = 8.0; Z = 0.0;}

    under Node[10] {X = 0.0; Y = 0.0; Z = 0.0;}
    under Node[11] {X = 0.0; Y = 2.0; Z = 0.0;}
    under Node[12] {X = 0.0; Y = 4.0; Z = 0.0;}
    under Node[13] {X = 0.0; Y = 6.0; Z = 0.0;}
    under Node[14] {X = 0.0; Y = 8.0; Z = 0.0;}
  
    under Node[15] {X = 2.0; Y = 0.0; Z = 0.0;}
    under Node[16] {X = 2.0; Y = 2.0; Z = 0.0;}
    under Node[17] {X = 2.0; Y = 4.0; Z = 0.0;}
    under Node[18] {X = 2.0; Y = 6.0; Z = 0.0;}
    under Node[19] {X = 2.0; Y = 8.0; Z = 0.0;}
  
    under Node[20] {X = 4.0; Y = 0.0; Z = 0.0;}
    under Node[21] {X = 4.0; Y = 2.0; Z = 0.0;}
    under Node[22] {X = 4.0; Y = 4.0; Z = 0.0;}
    under Node[23] {X = 4.0; Y = 6.0; Z = 0.0;}
    under Node[24] {X = 4.0; Y = 8.0; Z = 0.0;}

  } 
  VerboseScreenOutput = Low; //High;
}

