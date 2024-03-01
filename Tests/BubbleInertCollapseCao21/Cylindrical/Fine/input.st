// units: mm, g, s
// Ref: Cao et al., JFM 2021 (Figs. 6, 7)

under Mesh {
  Type = Cylindrical;
  X0   = -20.0;
  Xmax = 20.0;
  Y0   = 0.0;
  Ymax = 20.0;
  Z0   = -0.001;
  Zmax = 0.001;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate = 0.0;  CellWidth = 7.5e-3;}
  under ControlPointX[1] {Coordinate = 0.7;  CellWidth = 7.5e-3;}
  under ControlPointX[2] {Coordinate = 2.0;  CellWidth = 1.8e-2;}
  under ControlPointX[3] {Coordinate = 5.0;  CellWidth = 4e-2;}
  under ControlPointX[4] {Coordinate = 10.0;  CellWidth = 2e-1;}
  under ControlPointX[5] {Coordinate = 20.0;  CellWidth = 4e-1;}
  under ControlPointX[6] {Coordinate = -0.7;  CellWidth = 7.5e-3;}
  under ControlPointX[7] {Coordinate = -2.0;  CellWidth = 1.8e-2;}
  under ControlPointX[8] {Coordinate = -5.0;  CellWidth = 4e-2;}
  under ControlPointX[9] {Coordinate = -10.0;  CellWidth = 2e-1;}
  under ControlPointX[10] {Coordinate = -20.0;  CellWidth = 4e-1;}

  under ControlPointY[0] {Coordinate = 0.0;  CellWidth = 7.5e-3;}
  under ControlPointY[1] {Coordinate = 0.7;  CellWidth = 7.5e-3;}
  under ControlPointY[2] {Coordinate = 2.0;  CellWidth = 1.8e-2;}
  under ControlPointY[3] {Coordinate = 5.0;  CellWidth = 4e-2;}
  under ControlPointY[4] {Coordinate = 10.0;  CellWidth = 2e-1;}
  under ControlPointY[5] {Coordinate = 20.0;  CellWidth = 4e-1;}
 
  NumberOfCellsZ = 1;
}

under Equations {
  under Material[0] { // liquid (Muller et al. Computer & Fluid, 2009)
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.1;
      PressureConstant = 2.045e9;
    }
  }
  under Material[1] { // gas bubble  (air)
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
    } 
    PressureCutOff = 1.0e-12;
  }
}

under InitialCondition {
  under GeometricEntities {
    under Sphere[0] {
      Center_x = 0.0;
      Center_y = 0.0;
      Center_z = 0.0;
      Radius = 0.7469;
      under InitialState {
        MaterialID = 1;
        Density = 0.957e-9; 
        VelocityX = 0.0;
        VelocityY = 0.0;
        VelocityZ = 0.0;
        Pressure = 100.0;
      }
    }
  }
}

under BoundaryConditions {
  under Farfield {
    MaterialID = 0;
    Density = 0.0010; 
    VelocityX = 0.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
    Pressure = 1.0e5;
  }
}

under Space {
  under NavierStokes {
    Flux = LocalLaxFriedrichs; //HLLC; //Roe;
    under Reconstruction {
      Type = Linear;    
      VariableType = Primitive; //Conservative; //ConservativeCharacteristic;
      Limiter = GeneralizedMinMod;
      GeneralizedMinModCoefficient = 1.0; //default: 1.2 (larger --> less dissipation)

      under Fixes {
        under Sphere[0] {
          Center_x = 0.0;
          Center_y = 0.0;
          Center_z = 0.0;
          Radius = 0.1;
        }
      }
    }
  }

  under LevelSet[0] {
    MaterialID = 1;
    Solver = FiniteDifference;

    BoundaryConditionX0   = LinearExtrapolation;
    BoundaryConditionXmax = LinearExtrapolation;
    BoundaryConditionY0   = ZeroNeumann;
    BoundaryConditionYmax = ZeroNeumann;
    BoundaryConditionZ0   = ZeroNeumann;
    BoundaryConditionZmax = ZeroNeumann;
  }
}

under MultiPhase {
  Flux = Numerical;
  ReconstructionAtInterface = Constant;
  PhaseChange = RiemannSolution;
  RiemannNormal = Mesh;
}

under Time {
  Type = Explicit;
  MaxTime = 1.2e-4;
  CFL = 0.4;
  under Explicit {
    Type = RungeKutta2;
  }
}


under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 5.0e-6;
  Velocity = On;
  LevelSet0 = On;
  MaterialID = On;
  Pressure = On;
  Density = On;

  VerboseScreenOutput = Medium;

  under MaterialVolumes { //NOTE: does not account for spherical symmetry!
    FileName = "volumes";
    TimeInterval = 1.0e-7;
  }
}

under TerminalVisualization {
  Plane = XY;
  Variable = Pressure; //Velocity;
  ClockTimeInterval = 20;
  HorizontalMin = -2.0;
  HorizontalMax = 2.0;
  VerticalMin = -0.0;
  VerticalMax = 2.0;
}


