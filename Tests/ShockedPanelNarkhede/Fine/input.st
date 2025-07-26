// Units: g, mm, s

under ConcurrentPrograms {
  under AeroS {
    FSIAlgorithm = ByAeroS;
  }
}

under Mesh {
  Type = ThreeDimensional;
  X0   =  -50.0;
  Xmax =  250.0;
  Y0   =    0.0;
  Ymax =   65.0;
  Z0   =   -0.8;
  Zmax =    0.8;

  NumberOfCellsZ = 1;

  under ControlPointX[0] {Coordinate = -50.0; CellWidth  = 1.0;}
  under ControlPointX[1] {Coordinate = -20.0; CellWidth  = 0.5;}
  under ControlPointX[2] {Coordinate = -10.0; CellWidth  = 0.10;}
  under ControlPointX[3] {Coordinate =  10.0; CellWidth  = 0.10;}
  under ControlPointX[4] {Coordinate =  20.0; CellWidth  = 0.5;}
  under ControlPointX[5] {Coordinate = 250.0; CellWidth  = 1.0;}

  under ControlPointY[0] {Coordinate =  0.0; CellWidth  = 0.10;}
  under ControlPointY[1] {Coordinate = 65.0; CellWidth  = 0.10;}

  BoundaryConditionX0   = Inlet;
  BoundaryConditionXmax = SlipWall;
  BoundaryConditionY0   = SlipWall;
  BoundaryConditionYmax = SlipWall;
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

  under GeometricEntities {
    under Plane[0] {
      Point_x = -20.0;
      Normal_x = -1; // behind the plane
      
      under InitialState {
        MaterialID = 0;
        Density = 1.63e-6;
        VelocityX = 1.096e5;
        Pressure = 1.557e5;
      }

    }  
    under Plane[1] {
      Point_x = -20;
      Normal_x = 1; // front of the plane

      under InitialState {
        MaterialID = 0;
        Density = 1.2e-6;
        VelocityX = 0.0;
        Pressure = 1.01e5;   
      }
    }
  }

}

under BoundaryConditions {
  under Inlet {
    MaterialID = 0;
    Density = 1.63e-6;
    VelocityX = 1.096e5;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.557e5;
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
  MaxTime = 5.0e-1;
  CFL    = 0.4;
  under Explicit {
    Type = RungeKutta3; 
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  Frequency = 500;
  Density = On;
  Velocity = On;
  Pressure = On;
  InternalEnergyPerUnitMass = On;
  MaterialID = On;

  //MeshPartition = "partition";

  under Probes {
    //TimeInterval = 2.0e-6;
    Frequency = 50;
    Pressure = "pressure_probes.txt";
    VelocityX = "vx_probes.txt";
    VelocityY = "vy_probes.txt";
    VelocityZ = "vz_probes.txt";
    Density  = "density_probes.txt";
    under Node[0] {X = -10.0; Y =  65.0; Z =   0.0;}
  }

  VerboseScreenOutput = Medium; //Low //High
}
