//units: g, mm, s, K

//Overlaying meshes ("Chimera")
under ConcurrentPrograms {
  under M2CTwin {
    Type = OversetGrids;
  }
}

under Mesh {

  Type = Cylindrical;
  X0   = -8.0;
  Xmax = 8.0;
  Y0   = 0.0;
  Ymax = 8.0;
  Z0   = -0.00125;
  Zmax = 0.00125;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate = -8.0; CellWidth  = 0.3;}
  under ControlPointX[1] {Coordinate = -6.0; CellWidth  = 0.2;}
  under ControlPointX[2] {Coordinate = -4.0; CellWidth  = 0.10;}
  under ControlPointX[3] {Coordinate = -2.0; CellWidth  = 0.03;}
  under ControlPointX[4] {Coordinate =  2.0; CellWidth  = 0.03;}
  under ControlPointX[5] {Coordinate =  4.0; CellWidth  = 0.10;}
  under ControlPointX[6] {Coordinate =  6.0; CellWidth  = 0.2;}
  under ControlPointX[7] {Coordinate =  8.0; CellWidth  = 0.3;}

  under ControlPointY[0] {Coordinate =  0.0; CellWidth  = 0.03;}
  under ControlPointY[1] {Coordinate =  2.0; CellWidth  = 0.03;}
  under ControlPointY[2] {Coordinate =  4.0; CellWidth  = 0.10;}
  under ControlPointY[3] {Coordinate =  6.0; CellWidth  = 0.2;}
  under ControlPointY[4] {Coordinate =  8.0; CellWidth  = 0.3;}

  NumberOfCellsZ = 1;
}

under Equations {
  under Material[0] { //water (Coralic and Colonius, 2014)
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 6.12;
      PressureConstant = 3.43e8;

      ReferenceTemperature = 273.15;

      // use cp and h (the cp and cv of water are similar up to 100C)
      SpecificHeatAtConstantPressure = 4.2e9;  //mm2/(s2.K)
      ReferenceSpecificEnthalpy      = 4.1011e11; //h0 = e0 + p0/rho0, based on p0=1.0e5, rho0 = 1.0e-3
                                                  //Note: e0 (~4e11) >> p0/rho0 (1e8)

      DensityCutOff  = 1.0e-12;
      PressureCutOff = -3.0e8;
    }
  }
}


under BoundaryConditions {
  under Farfield {
    Density = 0.0010; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}


under Space {
  under NavierStokes {
    Flux = LocalLaxFriedrichs;
    under Reconstruction {
      Type = Linear;    
      VariableType = ConservativeCharacteristic;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.1;
    }
  }
}

under Time {
  Type = Explicit;
  MaxTime = 5.0e-5;
  CFL = 0.4;
  under Explicit {
    Type = RungeKutta2; 
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution_2";
  TimeInterval = 1.0e-7;
  Velocity = On;
  Pressure = On;
  Density = On;

  //MeshPartition = "mesh2";

  VerboseScreenOutput = Low;//High;
}
