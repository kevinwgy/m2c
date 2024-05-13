// A trivial case (uniform flow, non-dimensional)

under Mesh {

  Type = ThreeDimensional;
  X0   = 0.0;
  Xmax = 10.0;
  Y0   = 0.0;
  Ymax = 10.0;
  Z0   = -0.5;
  Zmax = 0.5;

  BoundaryConditionX0   = Inlet;
  BoundaryConditionXmax = Inlet2;
  BoundaryConditionY0   = Symmetry;
  BoundaryConditionYmax = Symmetry;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

/*
  under ControlPointX[1] {Coordinate = -4.0; CellWidth  = 0.10;}
  under ControlPointX[2] {Coordinate = -2.0; CellWidth  = 0.025;}
  under ControlPointX[3] {Coordinate = -1.0; CellWidth  = 0.0125;}
  under ControlPointX[4] {Coordinate = -0.5; CellWidth  = 0.0035;}
  under ControlPointX[5] {Coordinate =  0.3; CellWidth  = 0.0035;}
  under ControlPointX[6] {Coordinate =  1.0; CellWidth  = 0.0095;}
  under ControlPointX[7] {Coordinate =  2.0; CellWidth  = 0.025;}
  under ControlPointX[8] {Coordinate =  4.0; CellWidth  = 0.10;}
  under ControlPointX[9] {Coordinate =  6.0; CellWidth  = 0.2;}

  under ControlPointY[0] {Coordinate =  0.0; CellWidth  = 0.0035;}
  under ControlPointY[1] {Coordinate =  0.4; CellWidth  = 0.0035;}
  under ControlPointY[2] {Coordinate =  1.0; CellWidth  = 0.0125;}
  under ControlPointY[3] {Coordinate =  2.0; CellWidth  = 0.025;}
  under ControlPointY[4] {Coordinate =  4.0; CellWidth  = 0.10;}
  under ControlPointY[5] {Coordinate =  6.0; CellWidth  = 0.2;}
*/

  NumberOfCellsX = 100;
  NumberOfCellsY = 100;
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
      DynamicViscosity = 1.0; //Re = rho*L*u/mu
    }
  }
}

under InitialCondition {
  under DefaultInitialState {

  }
}

under BoundaryConditions {
  under Inlet {
    VelocityX = -1.0;
    VelocityY = 0.0;
    VelocityZ = 0.0;
  }
  under Inlet2 {
    VelocityX = -1.0;
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
  MaxTime = 1.0;
  TimeStep = 0.1;
  //CFL = 500;

  under SemiImplicit {
    Type = SIMPLEC; //or SIMPLER or SIMPLEC
    //E = 5.0;
    E = 1.0e12;
    AlphaP = 0.8; 

    ConvergenceTolerance = 1.0e-4;
    MaxIts = 50;

    under LinearSolverForVelocity {
      Type = FlexibleGMRes; //ImprovedStabilizedBiCG;
      Preconditioner = BlockJacobi;
      RelativeErrorTolerance = 1.0e-5;
      MaxIts = 500;
      //PrintLogToScreen = Yes;
    }

    under LinearSolverForPressure {
      Type = FlexibleGMRes; //ImprovedStabilizedBiCG;
      Preconditioner = BlockJacobi;
      RelativeErrorTolerance = 1.0e-5;
      MaxIts = 500;
      //PrintLogToScreen = Yes;
    }
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  Frequency = 1;
  //TimeInterval = 2.0e-8;
  Velocity = On;
  Pressure = On;
  VerboseScreenOutput = Medium;
}
