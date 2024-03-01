// units: mm, g, s, K, (Pa)

under Mesh {
  Type = Cylindrical;
  X0   = 0.0;
  Xmax = 1.5;
  Y0   = 0.0;
  Ymax = 1.5;
  Z0   = -1.0e-3;
  Zmax = 1.0e-3; 

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  NumberOfCellsX = 150;
  NumberOfCellsY = 150;
  NumberOfCellsZ = 1;
}

under Equations {
  under Material[0] {
    EquationOfState = StiffenedGas;
    under StiffenedGasModel { //water
      SpecificHeatRatio = 6.12;
      PressureConstant = 3.43e8;
    }
    under HyperelasticityModel {
      //Type = ModifiedSaintVenantKirchhoff;
      //Type = SaintVenantKirchhoff;
      Type = NeoHookean;
      //Type = MooneyRivlin;
      StressTensor = Full;
      YoungsModulus = 1.0e7; //Pa "soft tissue"
      PoissonsRatio = 0.495; //close to 0.5
      C01 = 8.361204013377926e+05; // 1/4 of shear modulus (C01 + C10 = mu/2)
    }
  }
}

under Time {
  MaxIts = 0;
}

under Output {
  Prefix = "results/";
  under Probes {
    under Node[0] {
      X = 7.550000e-01;
      Y = 2.050000e-01;
      Z = 0.0;
    }
    PrincipalElasticStresses = "probe_stresses_2Dcyl_NH.txt";
  }
}

