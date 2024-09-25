// Flow over a flat plate
// Re = 5e6
under Mesh {

  Type = ThreeDimensional; //Parameters from https://turbmodels.larc.nasa.gov/flatplate_sa.html
  X0   = -0.5;
  Xmax = 2.0;
  Y0   = 0.0;
  Ymax = 1.0;
  Z0   = -0.005;
  Zmax = 0.005;

  BoundaryConditionX0   = Inlet;
  BoundaryConditionXmax = Outlet;
  BoundaryConditionY0   = Symmetry; //Wall starts somewhere (specified below)
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate = -0.5;  CellWidth  = 5e-2;}
  under ControlPointX[1] {Coordinate = -0.3;  CellWidth  = 1e-2;}
  under ControlPointX[2] {Coordinate = -0.1;  CellWidth  = 1e-3;}
  under ControlPointX[3] {Coordinate = 0.1;  CellWidth  = 1e-3;}
  under ControlPointX[4] {Coordinate = 0.5;  CellWidth  = 1e-2;}
  under ControlPointX[5] {Coordinate = 1.0;  CellWidth  = 5e-2;}

  under ControlPointY[0] {Coordinate = 0.0;   CellWidth  = 1.0e-4;}
  under ControlPointY[1] {Coordinate = 0.01;   CellWidth  = 3.0e-4;}
  under ControlPointY[2] {Coordinate = 0.02;   CellWidth  = 5.0e-4;} 
  under ControlPointY[3] {Coordinate = 0.03;  CellWidth  = 6.0e-4;}
  under ControlPointY[4] {Coordinate = 0.05;  CellWidth  = 9.0e-4;}
  under ControlPointY[5] {Coordinate = 0.1;  CellWidth  = 5.0e-3;}
  under ControlPointY[6] {Coordinate = 0.2;   CellWidth  = 5.0e-2;}
  under ControlPointY[7] {Coordinate = 0.5;   CellWidth  = 5.0e-2;}
  //under ControlPointY[6] {Coordinate = 2.0;   CellWidth  = 5.0e-2;} //addition of control point at boundary may help according to M2C documentation

//  NumberOfCellsX = 100;
//  NumberOfCellsY = 100;
  NumberOfCellsZ = 1;
}


under Equations {
  under Material[0] {
    EquationOfState = HomogeneousIncompressible;
    under HomogeneousIncompressibleModel {
      Density = 1.0;
      ReferencePressure = 0.0;
    }
    under ViscosityModel {
      Type = Constant;
      DynamicViscosity = 2.0e-7; //Re = rho*L*u/mu
    }
  }
}

under RANSTurbulence {
  Model = SpalartAllmaras;
  ExampleProblem = FlatPlate;
  ExampleParameterReal1 = 0.0;
  NuTildeAtFarfield = 8.0e-7; 

  //CW1ReductionFactor = 10.0;
  //CW1ReductionTime1  = 0.4955597;
  //CW1ReductionTime2  = 2.495356;
}


under InitialCondition {
  under DefaultInitialState {
    VelocityX = 1.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
  }
}

under BoundaryConditions {
  under Inlet {
    VelocityX = 1.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
  }
}

under Space {
  under NavierStokes {
  }
}

under Time {
  Type = SemiImplicit;
  LocalTimeStepping = Off;

  //MaxTime = 0.5;
  CFL = 1.0; 
  ConvergenceTolerance = 2.0e-7;

  under SemiImplicit {
    Type = SIMPLEC;
    E = 5.0;
    AlphaP = 0.8; 
    FixPressureAtOneCorner = Yes;

    ConvergenceTolerance = 1.0e-5;
    MaxIts = 20;

    under LinearSolverForVelocity {
      Type = StabilizedBiCG;
      Preconditioner = BlockJacobi;
      RelativeErrorTolerance = 1.0e-7;
      MaxIts = 500;
      //PrintLogToScreen = Yes;
    }

    under LinearSolverForPressure {
      Type = StabilizedBiCG;
      Preconditioner = BlockJacobi;
      RelativeErrorTolerance = 1.0e-7;
      MaxIts = 2000;
      //PrintLogToScreen = Yes;
      //LogFile = "pressure_lin.out";
    }

    under LinearSolverForTurbulenceClosure {
      Type = FlexibleGMRes;
      Preconditioner = BlockJacobi;
      RelativeErrorTolerance = 1.0e-7;
      MaxIts = 700;
      //PrintLogToScreen = Yes;
    }
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  Frequency = 5;
  //TimeInterval = 0.005;
  Density = On;
  Velocity = On;
  Pressure = On;
  KinematicEddyViscosity = On;
  VerboseScreenOutput = Medium;

  /*under LinePlot[0] {
    FileName = "vertical_x0p97";
    NumberOfPoints = 10000;
    Frequency = 800;
    X0   = 0.97;
    Y0   = 0.0; 
    Z0   = 0.0;
    Xmax = 0.97;
    Ymax = 0.5;
    Zmax = 0.0;
  }*/  
}

/*under TerminalVisualization {
  Plane = XY;
  Variable = Velocity;
  ClockTimeInterval = 10;
}*/
