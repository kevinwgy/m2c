// units: mm, g, s, K, (Pa)

under Mesh {
  X0   = 0.0;
  Xmax = 0.5;
  Y0   = 0.0;
  Ymax = 0.5;
  Z0   = 0.0;
  Zmax = 1.5;

  BoundaryConditionX0   = Symmetry;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Farfield;
  BoundaryConditionZmax = Farfield;

  NumberOfCellsX = 50;
  NumberOfCellsY = 50;
  NumberOfCellsZ = 150;
}

under Equations {
  under Material[0] {
    EquationOfState = StiffenedGas;
    under StiffenedGasModel { //water
      SpecificHeatRatio = 6.12;
      PressureConstant = 3.43e8;
    }
    under HyperelasticityModel {
      //Type = SaintVenantKirchhoff;
      //Type = ModifiedSaintVenantKirchhoff;
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
      X = 2.050000e-01;
      Y = 5.000000e-03;
      Z = 7.550000e-01;
    }
    under Node[1] {
      X = 1.450000e-01;
      Y = 1.450000e-01;
      Z = 7.550000e-01;
    }
    under Node[2] {
      X = 5.000000e-03;
      Y = 2.050000e-01;
      Z = 7.550000e-01;
    }
    PrincipalElasticStresses = "probe_stresses_3d_NH.txt";
  }
}

