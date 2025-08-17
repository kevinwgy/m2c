//Units: mm, g, s

under Mesh {
  Type = Cylindrical;
  X0   = -8.0;
  Xmax = 16.0;
  Y0   = 0.0;
  Ymax = 6.5;
  Z0   = -0.015;
  Zmax = 0.015;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate =   -8.0; CellWidth  = 0.1;}
  under ControlPointX[1] {Coordinate =   -3.0; CellWidth  = 0.008;}
  under ControlPointX[2] {Coordinate =    10.0; CellWidth  = 0.008;}
  under ControlPointX[3] {Coordinate =   16.0; CellWidth  = 0.1;}

  under ControlPointY[0] {Coordinate =    0.0; CellWidth  = 0.008;}
  under ControlPointY[1] {Coordinate =    3.0; CellWidth  = 0.008;}
  under ControlPointY[2] {Coordinate =    6.5; CellWidth  = 0.1;}

  NumberOfCellsZ = 1;
}

under Equations {
  under Material[0] { //Air
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0;

      // Apply a simple temperature law, not "thermodynamically consistent"
      SpecificHeatAtConstantVolume = 717.1e6;  // 717.1 J/(K.kg) = 717.1e6 mm2.s2/K
      ReferenceTemperature = 300.0; //K
      ReferenceSpecificInternalEnergy = 2.083333e11; //1.0e5/((1.4-1.0)*1.2e-6) ensures T = T0 in ref. cond.
    }
    //DensityCutOff = 1.0e-12;
    //PressureCutOff = 1.0;
    //DensityPrescribedAtFailure = 1.2e-6; //density of air
  }

  under Material[1] { //Helium
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.66;
      PressureConstant = 0;

      // Apply a simple temperature law, not "thermodynamically consistent" (WARNING: THESE NUMBERS ARE WRONG)
      SpecificHeatAtConstantVolume = 717.1e6;  // 717.1 J/(K.kg) = 717.1e6 mm2.s2/K
      ReferenceTemperature = 300.0; //K
      ReferenceSpecificInternalEnergy = 2.083333e11; //1.0e5/((1.4-1.0)*1.2e-6) ensures T = T0 in ref. cond.
    }
    //DensityCutOff = 1.0e-12;
    //PressureCutOff = 1.0;
    //DensityPrescribedAtFailure = 1.2e-6; //density of air
  }
}


under InitialCondition {

  under GeometricEntities {
    under Point[0] {
      X = 3.0;
      Y = 0.2;
      Z = 0.0;
      under InitialState {
        MaterialID = 1; //helium
        Density = 2.5e-6;
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 4.0e5;
      }
    }

    under Point[1] {
      X = -3.5;
      Y = 0.2;
      Z = 0.0;
      under InitialState {
        MaterialID = 0; //air
        Density = 2.0e-6;
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 4.0e5;
      }
    }
  }

}

under BoundaryConditions {
  under Farfield {
    Density = 1.2e-6; //density of air
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}

under EmbeddedBoundaryMethod {
  RiemannNormal = Average; // Mesh; //EmbeddedSurface;
  ReconstructionAtInterface = Constant;
  under EmbeddedSurfaces {
    under Surface[0] {
      MeshFile = "Geometry/box2d.top";
      SurfaceThickness = 1.0e-5;
      GaussQuadrature = OnePoint;
      GaussPointsLofting = 0.0;
      InternalPressure = 1.0e5;
      //ConstantReconstructionDepth = 0.0;

      under Output {
        TimeInterval = 2.0e-7;
        Prefix = "results/";
        Mesh = "target.top";
        Displacement = "target_disp.txt";
        Solution = "target_load.txt";
        ContactSurfaceOutput = "target_contact.top";
      }
    }
    under Surface[1] {
      MeshFile = "Geometry/sphere2d_out.top";
      SurfaceThickness = 1.0e-5;
      GaussQuadrature = OnePoint;
      GaussPointsLofting = 0.0;
      InternalPressure = 1.0e5;
      //ConstantReconstructionDepth = 0.0;

      UserDefinedDynamicsCalculator = "./UserDefinedDynamics.so";

      under Output {
        TimeInterval = 2.0e-7;
        Prefix = "results/";
        Mesh = "impactor.top";
        Displacement = "impactor_disp.txt";
        Solution = "impactor_load.txt";
        ContactSurfaceOutput = "impactor_contact.top";
      }
    }

    under SurfaceIntersection[0] {
      Surface1 = 0;
      Surface2 = 1;
      Enclosure = IgnoreSurface1;
    }
  }
}

under Space {
  under NavierStokes {
    Flux = LocalLaxFriedrichs; //HLLC;
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
  MaxTime = 2.9e-5;
  CFL    = 0.3;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 2.0e-7;
  Density = On;
  Velocity = On;
  Pressure = On;
  //Temperature = On;
  MaterialID = On;

  MeshInformation = "mesh.txt";

  VerboseScreenOutput = Medium; //Low; //High;
}

under TerminalVisualization {
  Plane = XY;
  Variable = MaterialID; //Pressure; //Velocity;
  ClockTimeInterval = 60;
  HorizontalMin = -2.0;
  HorizontalMax = 2.0;
  VerticalMin = 0.0;
  VerticalMax = 3.0;
}
