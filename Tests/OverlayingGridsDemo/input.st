//units: g, mm, s, K

//Overlaying meshes ("Chimera")
under ConcurrentPrograms {
  under M2CTwin {
    Type = OversetGrids;
  }
}

under Mesh {

  Type = Cylindrical;
  X0   = -1.5;
  Xmax = 1.5;
  Y0   = 0.0;
  Ymax = 1.5;
  Z0   = -0.00125;
  Zmax = 0.00125;

  BoundaryConditionX0   = Overset;
  BoundaryConditionXmax = Overset;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Overset;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate = -1.5; CellWidth  = 0.020;}
  under ControlPointX[1] {Coordinate = -1.0; CellWidth  = 0.0075;}
  under ControlPointX[2] {Coordinate = -0.5; CellWidth  = 0.0025;}
  under ControlPointX[3] {Coordinate =  0.3; CellWidth  = 0.0025;}
  under ControlPointX[4] {Coordinate =  1.0; CellWidth  = 0.0075;}
  under ControlPointX[5] {Coordinate =  1.5; CellWidth  = 0.020;}

  under ControlPointY[0] {Coordinate =  0.0; CellWidth  = 0.0025;}
  under ControlPointY[1] {Coordinate =  0.4; CellWidth  = 0.0025;}
  under ControlPointY[2] {Coordinate =  1.0; CellWidth  = 0.0075;}
  under ControlPointY[3] {Coordinate =  1.5; CellWidth  = 0.020;}

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

  under Material[1] { //water vapor
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.34;
      PressureConstant = 0.0;

      ReferenceTemperature = 373.15; 

      // use cp and h
      SpecificHeatAtConstantPressure = 2.0e9;  //1bar, 130 C, https://webbook.nist.gov/chemistry/fluid/
      ReferenceSpecificEnthalpy      = 3.0881e12; //= 4.1011e11 + 4.2e9*(373.15-273.15) + 2.258e12 (lat.heat)

      DensityCutOff  = 1.0e-13;
      PressureCutOff = 1.0e-13;
    }
  }

}



under EmbeddedBoundaryMethod {
  RiemannNormal = Mesh; //Average;
  ReconstructionAtInterface = Constant;
  under EmbeddedSurfaces {
    under Surface[0] {
      MeshFile = "laser_probe.top";
      SurfaceThickness = 1.0e-8;
      GaussQuadrature = OnePoint;
      GaussPointsLofting = 0.0;
      InternalPressure = 0.0;
      ConstantReconstructionDepth = 0.003;
      under Output {
        //Frequency = 1;
        TimeInterval = 1.0e-7;
        Prefix = "results/";
        Mesh = "surface.top";
        Displacement = "surface_disp.txt";
        Solution = "surface_load.txt";
      }
    }
  }
}


under InitialCondition {

  under DefaultInitialState {
    Density = 0.0010; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }

  under GeometricEntities {
    under Sphere[0] { //expanding bubble
      Center_x = -0.2;
      Center_y = 0.0;
      Center_z = 0;
      Radius = 0.1;
      under InitialState {
        MaterialID = 1;
        Density = 5.0e-5;
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 5.0e7;
      }
    }

  }
}

/*
under BoundaryConditions {
  under Farfield {
    Density = 0.0010; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}
*/

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


  under LevelSet[0] {
    MaterialID = 1;
    Solver = FiniteDifference;

    Bandwidth = 6;

    BoundaryConditionX0   = LinearExtrapolation;
    BoundaryConditionXmax = LinearExtrapolation;
    BoundaryConditionY0   = ZeroNeumann;
    BoundaryConditionYmax = LinearExtrapolation;
    BoundaryConditionZ0   = ZeroNeumann;
    BoundaryConditionZmax = ZeroNeumann;

    under Reinitialization {
      Frequency = 2;
      MaxIts = 400;
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
  MaxTime = 2.0e-6;
  CFL = 0.30;
  under Explicit {
    Type = RungeKutta2; 
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 1.0e-7;
  Velocity = On;
  LevelSet0 = On;
  MaterialID = On;
  Pressure = On;
  Density = On;

  VerboseScreenOutput = Medium; //Low;//High;
}
