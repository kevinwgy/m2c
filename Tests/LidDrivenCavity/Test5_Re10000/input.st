// Lid-driven cavity
// Re = 10000
// https://www.comsol.com/blogs/how-to-solve-a-classic-cfd-benchmark-the-lid-driven-cavity-problem/

under Mesh {

  Type = ThreeDimensional;
  X0   = 0.0;
  Xmax = 1.0;
  Y0   = 0.0;
  Ymax = 1.0;
  Z0   = -0.005;
  Zmax = 0.005;

  BoundaryConditionX0   = StickWall;
  BoundaryConditionXmax = StickWall;
  BoundaryConditionY0   = StickWall;
  BoundaryConditionYmax = Inlet;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate = 0.0;  CellWidth  = 0.001;}
  under ControlPointX[1] {Coordinate = 0.2; CellWidth  = 0.0025;}
  under ControlPointX[2] {Coordinate = 0.5; CellWidth  = 0.005;}
  under ControlPointX[3] {Coordinate = 0.8; CellWidth  = 0.0025;}
  under ControlPointX[4] {Coordinate = 1.0;  CellWidth  = 0.001;}

  under ControlPointY[0] {Coordinate = 0.0;  CellWidth  = 0.001;}
  under ControlPointY[1] {Coordinate = 0.2; CellWidth  = 0.0025;}
  under ControlPointY[2] {Coordinate = 0.5; CellWidth  = 0.005;}
  under ControlPointY[3] {Coordinate = 0.8; CellWidth  = 0.0025;}
  under ControlPointY[4] {Coordinate = 1.0;  CellWidth  = 0.001;}

//  NumberOfCellsX = 100;
//  NumberOfCellsY = 100;
  NumberOfCellsZ = 1;
}


under Equations {
  under Material[0] {
    EquationOfState = HomogeneousIncompressible;
    under HomogeneousIncompressibleModel {
      Density = 1.0;
      ReferencePressure = 1.0;
    }
    under ViscosityModel {
      Type = Constant;
      DynamicViscosity = 0.0001; //Re = rho*L*u/mu
    }
  }
}

under InitialCondition {
  under DefaultInitialState {
    VelocityX = 0.0;
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
  LocalTimeStepping = On;
  //MaxTime = 0.5;
  CFL = 1000.0; 
  ConvergenceTolerance = 2.0e-7;

  under SemiImplicit {
    Type = SIMPLEC; //or SIMPLER or SIMPLEC or PISO
    E = 5.0;
    //E = 1.0e12;
    AlphaP = 0.8; 

    ConvergenceTolerance = 1.0e-7;
    MaxIts = 50;

    under LinearSolverForVelocity {
      //Type = FlexibleGMRes; //ImprovedStabilizedBiCG;
      Type = ImprovedStabilizedBiCG;
      Preconditioner = BlockJacobi;
      RelativeErrorTolerance = 1.0e-7;
      MaxIts = 500;
      //PrintLogToScreen = Yes;
    }

    under LinearSolverForPressure {
      //Type = FlexibleGMRes; //ImprovedStabilizedBiCG;
      Type = ImprovedStabilizedBiCG;
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
  Frequency = 200;
  //TimeInterval = 0.005;
  Density = On;
  Velocity = On;
  Pressure = On;
  VerboseScreenOutput = Medium;

  under LinePlot[0] {
    FileName = "vertical_centerline";
    NumberOfPoints = 300;
    Frequency = 800;
    X0   = 0.5;
    Y0   = 0.0; 
    Z0   = 0.0;
    Xmax = 0.5;
    Ymax = 1.0;
    Zmax = 0.0;
  }  

  under LinePlot[1] {
    FileName = "horizontal_centerline";
    NumberOfPoints = 300;
    Frequency = 800;
    X0   = 0.0; 
    Y0   = 0.5;
    Z0   = 0.0;
    Xmax = 1.0;
    Ymax = 0.5;
    Zmax = 0.0;
  }  
}

under TerminalVisualization {
  Plane = XY;
  Variable = Velocity;
  ClockTimeInterval = 10;
}
