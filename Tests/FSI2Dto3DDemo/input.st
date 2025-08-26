// Units: g, mm, s

under ConcurrentPrograms {
  under AeroS {
    FSIAlgorithm = ByAeroS;
  }
}

under Mesh {
  Type = Cylindrical;
  X0   = 0;
  Xmax = 1200;
  Y0   = 0;
  Ymax = 1200;
  Z0   = -2.5;
  Zmax =  2.5;

  NumberOfCellsZ = 1;

  under ControlPointX[0] {Coordinate =   0.0; CellWidth  = 5.0;}
  under ControlPointX[1] {Coordinate = 600.0; CellWidth  = 5.0;}
  under ControlPointX[2] {Coordinate = 1200.0; CellWidth  = 40.0;}

  under ControlPointY[0] {Coordinate =   0.0; CellWidth  = 5.0;}
  under ControlPointY[1] {Coordinate = 600.0; CellWidth  = 5.0;}
  under ControlPointY[2] {Coordinate = 1200.0; CellWidth  = 40.0;}

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
        Frequency = 1;
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

  //...Material 0: air
  under Material[0]{
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
  
  ApplyUserDataBeforeGeometricEntities = Yes; //default is No.

  UserDataFile = "PlanarWave.txt";

  under GeometricEntities {

    under Point[0] { //inside structure 
      X = 10.0;
      Y = 10.0;
      Z = 0.0;
      under InitialState {
        MaterialID = 0;
        Density = 1.177e-6;
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
      GeneralizedMinModCoefficient = 1.2;
    }
  }

}


under Time {
  Type = Explicit;
  MaxTime = 1.0e-2;
  CFL    = 0.4;
  under Explicit {
    Type = RungeKutta3; 
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  Frequency = 1;
  //TimeInterval = 1.0e-5;
  Density = On;
  Velocity = On;
  Pressure = On;
  MaterialID = On;

  //MeshPartition = "partition";

/*
  under Probes {
    //TimeInterval = 2.0e-6;
    Frequency = 1;
    Pressure = "pressure_probes.txt";
    VelocityX = "vx_probes.txt";
    VelocityY = "vy_probes.txt";
    VelocityZ = "vz_probes.txt";
    Density  = "density_probes.txt";
    under Node[0] {X = -485.0; Y = 0.0; Z = 0.0;}
    under Node[1] {X = 0.0; Y = 0.0; Z = -315.0;}
    under Node[2] {X = -485.0; Y = 0.0; Z = -315.0;}
  }
*/

  VerboseScreenOutput = Medium; //Low //High
}
