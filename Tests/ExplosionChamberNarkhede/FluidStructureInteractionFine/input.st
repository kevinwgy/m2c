// Units: g, mm, s

under ConcurrentPrograms {
  under AeroS {
    FSIAlgorithm = ByAeroS;
  }
}

under Mesh {
  Type = Cylindrical;
  X0   =    0;
  Xmax = 1200;
  Y0   =    0;
  Ymax = 1200;
  Z0   = -3.0;
  Zmax =  3.0;

  NumberOfCellsZ = 1;

  under ControlPointX[0] {Coordinate =    0.0; CellWidth  =  0.5;}
  under ControlPointX[1] {Coordinate =  300.0; CellWidth  =  0.5;}
  under ControlPointX[2] {Coordinate =  600.0; CellWidth  =  5.0;}
  under ControlPointX[3] {Coordinate = 1200.0; CellWidth  = 50.0;}

  under ControlPointY[0] {Coordinate =    0.0; CellWidth  =  0.5;}
  under ControlPointY[1] {Coordinate =  300.0; CellWidth  =  0.5;}
  under ControlPointY[2] {Coordinate =  600.0; CellWidth  =  5.0;}
  under ControlPointY[3] {Coordinate = 1200.0; CellWidth  = 50.0;}

  BoundaryConditionX0   = Symmetry;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;
}


under EmbeddedBoundaryMethod {
  RiemannNormal = Mesh; //Average;
  ReconstructionAtInterface = Constant;
  under EmbeddedSurfaces {
    under Surface[0] {
      SurfaceProvidedByAnotherSolver = Yes;
      SurfaceThickness = 1.0e-6;
      GaussQuadrature = OnePoint;
      GaussPointsLofting = 0.2; //non-D, w.r.t. local element size
      InternalPressure = 0.0;
      ConstantReconstructionDepth = 6.0; //dimensional (distance)
      //TwoDimensionalToThreeDimensionalMapping = NearestNeighbor;
      under Output {
        Frequency = 40;
        Prefix = "results/";
        Mesh = "surface.top";
        Displacement = "surface_disp.txt";
        Solution = "surface_load.txt";
        ContactSurfaceOutput = "surface_contact.top";
      }
    }
  }
}


under Equations {

  //...Material 0: burnt gas and ambient air
  under Material[0] {
    EquationOfState = JonesWilkinsLee;
    under JonesWilkinsLeeModel {
      Omega = 0.28;
      A1    = 3.712e11;
      A2    = 3.21e9;
      R1    = 4.15;
      R2    = 0.95;
      Rho0  = 1.63e-3;
    }
    PressureCutOff = 1.0e-2;
    DensityCutOff = 1.0e-11;
  }

  //...Material 1: external air
  under Material[1] {
    EquationOfState = StiffenedGas;
    under StiffenedGasModel{
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
    }
    PressureCutOff = 1.0;
    DensityCutOff  = 1.0e-11;
  }
}

under InitialCondition {
  
  // AN: Use this when starting from burning.
  //UserDefinedStateCalculator = "./StateCalculator/StateCalculator.so";

  UserDataFile = "SphericalShock.txt";

  under GeometricEntities {

      // Point inside the structure.
      under Point[0] {
        X = 0.0;
        Y = 10.0;
        Z = 0.0;

      // initialize with ambient air.
      under InitialState {
        MaterialID = 0;
        Density = 1.177e-6; 
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 1e5;
      }
    }

  }
}

under BoundaryConditions {
  under Farfield {
    MaterialID = 1;
    Density = 1.177e-6;
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5; 
  }
}

under Space {

  under NavierStokes {
    Flux = LocalLaxFriedrichs; //HLLC; //LocalLaxFriedrichs; //Roe; 
    under Reconstruction {
      Type = Linear;    
      VariableType = ConservativeCharacteristic;
      SlopeNearInterface = Zero;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.0;
    }
  }

}

under Time {
  Type = Explicit;
  MaxTime = 1.0e-2;
  CFL    = 0.4;
  under Explicit {
    Type = RungeKutta2; 
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  Frequency = 400;
  //TimeInterval = 2e-5;
  Density = On;
  Velocity = On;
  Pressure = On;
  MaterialID = On;
  InternalEnergyPerUnitMass = On;
  //MeshPartition = "partition";

  under EnergyIntegration {
    Frequency = 50;
    //TimeInterval = 1e-8;
    TotalEnergy = "total_energy.txt";
    KineticEnergy = "kinetic_energy.txt";
    InternalEnergy = "internal_energy.txt";
  }

  under Probes {
    //TimeInterval = 1e-8;
    Frequency = 1;
    Pressure = "pressure_probes.txt";
    VelocityX = "vx_probes.txt";
    VelocityY = "vy_probes.txt";
    VelocityZ = "vz_probes.txt";
    Density  = "density_probes.txt";
    MaterialID = "material_probes.txt";
// 1 mm offset sensors
    under Node[0] {X = 179.000; Y =   1.000; Z =   0.0;}
    under Node[1] {X = 161.049; Y =  75.661; Z =   0.0;}
    under Node[2] {X = 108.156; Y = 135.370; Z =   0.0;}
    under Node[3] {X =  30.000; Y = 159.000; Z =   0.0;}
    under Node[4] {X =  15.000; Y = 159.000; Z =   0.0;}
    under Node[5] {X =   1.000; Y = 159.000; Z =   0.0;}

  }
  VerboseScreenOutput = Medium; //Low //High
}

/*
under TerminalVisualization {
  Plane = XY;
  Variable = MaterialID;
  ClockTimeInterval = 10;
  HorizontalMin = 0;
  HorizontalMax = 350;
  VerticalMin = 0;
  VerticalMax = 350;
}
*/
