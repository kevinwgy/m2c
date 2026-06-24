// Units: g, mm, s
under Mesh {
  Type = ThreeDimensional;
  X0   = -300.0;
  Xmax =  300.0;
  Y0   = -300.0;
  Ymax =  300.0;
  Z0   =   -1.0;
  Zmax =    1.0;

  BoundaryConditionX0   = Farfield;
  BoundaryConditionXmax = Farfield;
  BoundaryConditionY0   = Farfield;
  BoundaryConditionYmax = Farfield;
  BoundaryConditionZ0   = Symmetry;
  BoundaryConditionZmax = Symmetry;

  under ControlPointX[0] {Coordinate =-300.0; CellWidth  = 20.0;}
  under ControlPointX[1] {Coordinate =-160.0; CellWidth  = 3.0;}
  under ControlPointX[2] {Coordinate = -90.0; CellWidth  = 0.7;}
  under ControlPointX[3] {Coordinate = -60.0; CellWidth  = 0.3;}
  under ControlPointX[4] {Coordinate =  60.0; CellWidth  = 0.3;}
  under ControlPointX[5] {Coordinate =  90.0; CellWidth  = 0.7;}
  under ControlPointX[6] {Coordinate = 160.0; CellWidth  = 3.0;}
  under ControlPointX[7] {Coordinate = 300.0; CellWidth  = 20.0;}

  under ControlPointY[0] {Coordinate =-300.0; CellWidth  = 20.0;}
  under ControlPointY[1] {Coordinate =-150.0; CellWidth  = 3.0;}
  under ControlPointY[2] {Coordinate = -90.0; CellWidth  = 0.7;}
  under ControlPointY[3] {Coordinate = -60.0; CellWidth  = 0.3;}
  under ControlPointY[4] {Coordinate =  60.0; CellWidth  = 0.3;}
  under ControlPointY[5] {Coordinate =  90.0; CellWidth  = 0.7;}
  under ControlPointY[6] {Coordinate = 150.0; CellWidth  = 3.0;}
  under ControlPointY[7] {Coordinate = 300.0; CellWidth  = 20.0;}

  NumberOfCellsZ = 1;
}

under Equations {
  under Material[0] { // air
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.4;
      PressureConstant = 0.0;
    }
  }
}

under InitialCondition {
  under GeometricEntities {
    under Sphere[0] {
      Center_x = 40.0;
      Center_y = 40.0;
      Center_z = 0.0;
      Radius = 20.0;
      under InitialState {
        MaterialID = 0;
        Density = 5.0e-6;
        VelocityX = -7.071e4;
        VelocityY = -7.071e4;
        VelocityZ = 0.0;
        Pressure = 1.200000e+06;
      }
    }
    under Sphere[1] {
      Center_x = -40.0;
      Center_y = 40.0;
      Center_z = 0.0;
      Radius = 20.0;
      under InitialState {
        MaterialID = 0;
        Density = 5.0e-6;
        VelocityX = 7.071e4;
        VelocityY = -7.071e4;
        VelocityZ = 0.0;
        Pressure = 1.120000e+06;
      }
    }
    under Sphere[2] {
      Center_x = -40.0;
      Center_y = -40.0;
      Center_z = 0.0;
      Radius = 20.0;
      under InitialState {
        MaterialID = 0;
        Density = 5.0e-6;
        VelocityX = 7.071e4;
        VelocityY = 7.071e4;
        VelocityZ = 0.0;
        Pressure = 1.200000e+06;
      }
    }
    under Sphere[3] {
      Center_x = 40.0;
      Center_y = -40.0;
      Center_z = 0.0;
      Radius = 20.0;
      under InitialState {
        MaterialID = 0;
        Density = 5.0e-6;
        VelocityX = -7.071e4;
        VelocityY = 7.071e4;
        VelocityZ = 0.0;
        Pressure = 1.200000e+06;
      }
    }
  }
}

under BoundaryConditions {
  under Farfield {
    Density = 1.225e-6;
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
      GeneralizedMinModCoefficient = 1.2;
    }
  }
}

under Time {
  Type = Explicit;
  MaxTime = 2.0e-4;
  CFL = 0.4;
  under Explicit {
    Type = RungeKutta2;
  }
}

under Output {
  Prefix = "results/";
  Solution = "solution";
  TimeInterval = 2.0e-5;
  Velocity = On;
  Pressure = On;
  Density = On;

  under Integration[0] {
    TimeInterval = 2.0e-6;    
    Volume      = "volume_0.txt";
    Mass        = "mass_0.txt";
    Momentum    = "momentum_0.txt";
    TotalEnergy = "totalenergy_0.txt";
    under Plane[0] {
      Point_x = 0.00001; //to avoid coinciding with a grid point
      Point_y = 0.00001;
      Point_z = 0.0;
      Normal_x = 1.0;
      Normal_y = 0.0;
      Normal_z = 0.0;
      Inclusion = Intersection;
      OperationOrder = 0; //!< 0 -> highest rank
    }
    under Plane[1]{
      Point_x = 0.00001; //to avoid coinciding with a grid point
      Point_y = 0.00001;
      Point_z = 0.0;
      Normal_x = 0.0;
      Normal_y = 1.0;
      Normal_z = 0.0;
      Inclusion = Intersection;
      OperationOrder = 1;
    }
  }
  under Integration[1] {
    TimeInterval = 2.0e-6;    
    Volume      = "volume_1.txt";
    Mass        = "mass_1.txt";
    Momentum    = "momentum_1.txt";
    TotalEnergy = "totalenergy_1.txt";
    under Plane[0] {
      Point_x = 0.00001; //to avoid coinciding with a grid point
      Point_y = 0.00001;
      Point_z = 0.0;
      Normal_x = -1.0;
      Normal_y = 0.0;
      Normal_z = 0.0;
      Inclusion = Intersection;
      OperationOrder = 0; //!< 0 -> highest rank
    }
    under Plane[1]{
      Point_x = 0.00001; //to avoid coinciding with a grid point
      Point_y = 0.00001;
      Point_z = 0.0;
      Normal_x = 0.0;
      Normal_y = 1.0;
      Normal_z = 0.0;
      Inclusion = Intersection;
      OperationOrder = 1;
    }
  }
  under Integration[2] {
    TimeInterval = 2.0e-6;    
    Volume      = "volume_2.txt";
    Mass        = "mass_2.txt";
    Momentum    = "momentum_2.txt";
    TotalEnergy = "totalenergy_2.txt";
    under Plane[0] {
      Point_x = 0.00001; //to avoid coinciding with a grid point
      Point_y = 0.00001;
      Point_z = 0.0;
      Normal_x = -1.0;
      Normal_y = 0.0;
      Normal_z = 0.0;
      Inclusion = Intersection;
      OperationOrder = 0; //!< 0 -> highest rank
    }
    under Plane[1]{
      Point_x = 0.00001; //to avoid coinciding with a grid point
      Point_y = 0.00001;
      Point_z = 0.0;
      Normal_x = 0.0;
      Normal_y = -1.0;
      Normal_z = 0.0;
      Inclusion = Intersection;
      OperationOrder = 1;
    }
  }
  under Integration[3] {
    TimeInterval = 2.0e-6;    
    Volume      = "volume_3.txt";
    Mass        = "mass_3.txt";
    Momentum    = "momentum_3.txt";
    TotalEnergy = "totalenergy_3.txt";
    under Plane[0] {
      Point_x = 0.00001; //to avoid coinciding with a grid point
      Point_y = 0.00001;
      Point_z = 0.0;
      Normal_x = 1.0;
      Normal_y = 0.0;
      Normal_z = 0.0;
      Inclusion = Intersection;
      OperationOrder = 0; //!< 0 -> highest rank
    }
    under Plane[1]{
      Point_x = 0.00001; //to avoid coinciding with a grid point
      Point_y = 0.00001;
      Point_z = 0.0;
      Normal_x = 0.0;
      Normal_y = -1.0;
      Normal_z = 0.0;
      Inclusion = Intersection;
      OperationOrder = 1;
    }
  }
}
