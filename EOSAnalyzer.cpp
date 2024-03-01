/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<EOSAnalyzer.h>
#include<fstream>
#include<cstring> //strcmp
#include<cassert>
using std::vector;

//------------------------------------------------------------------------------

EOSAnalyzer::EOSAnalyzer(ObjectMap<EOSTabulationData> &eos_tabulationMap_,
                         vector<VarFcnBase*>& vf_)
           : eos_tabulationMap(eos_tabulationMap_), vf(vf_) 
{ }

//------------------------------------------------------------------------------

EOSAnalyzer::~EOSAnalyzer()
{ }

//------------------------------------------------------------------------------

void
EOSAnalyzer::GenerateAllEOSTables()
{
  for(auto it = eos_tabulationMap.dataMap.begin(); 
           it != eos_tabulationMap.dataMap.end(); it++)
    GenerateEOSTable(it->first, *it->second);
}

//------------------------------------------------------------------------------

void
EOSAnalyzer::GenerateEOSTable(int table_id, EOSTabulationData& tab)
{

  int id = tab.materialid;
  if(id<0 || id>=(int)vf.size()) {
    print_error("*** Error: Cannot generate EOS table %d for undefined "
                "material id: %d.\n", table_id, id);
    exit(-1);
  }

  if(!strcmp(tab.filename,"")) {
    print_error("*** Error: Cannot generate EOS table %d for"
                   "material (%d). Output file is not specified.\n", table_id, id);
    exit(-1);
  }



  // Get Table Info.

  vector<vector<double> > result;
  double x0 = tab.x0, y0 = tab.y0, xmax = tab.xmax, ymax = tab.ymax;
  int Nx = tab.Nx, Ny = tab.Ny;
  string xfield, yfield, zfield;

  if(tab.output == EOSTabulationData::PRESSURE) {
    zfield = "Pressure";

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      xfield = "Density";
      yfield = "Specific Internal Energy";
      auto fun = [&](double rho, double e) {return vf[id]->GetPressure(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      xfield = "Density";
      yfield = "Temperature";
      auto fun = [&](double rho, double T) {
        double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T); 
        return vf[id]->GetPressure(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
    zfield = "Specific Internal Energy";

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      xfield = "Density";
      yfield = "Pressure";
      auto fun = [&](double rho, double p) {return vf[id]->GetInternalEnergyPerUnitMass(rho,p);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      xfield = "Density";
      yfield = "Temperature";
      auto fun = [&](double rho, double T) {return vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::DENSITY) {
    zfield = "Density";

    if(tab.xvar == EOSTabulationData::PRESSURE && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      xfield = "Density";
      yfield = "Specific Internal Energy";
      auto fun = [&](double p, double e) {return vf[id]->GetDensity(p,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::DP_DE) {
    zfield = "PressureDerivativeEnergy";

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      xfield = "Density";
      yfield = "Specific Internal Energy";
      auto fun = [&](double rho, double e) {return rho*vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      xfield = "Density";
      yfield = "Pressure";
      auto fun = [&](double rho, double p) {
        double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
        return rho*vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      xfield = "Density";
      yfield = "Temperature";
      auto fun = [&](double rho, double T) {
        double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T);
        return rho*vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::GRUNEISEN_PARAMETER) {
    zfield = "Gruneisen Parameter";

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      xfield = "Density";
      yfield = "Specific Internal Energy";
      auto fun = [&](double rho, double e) {return vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      xfield = "Density";
      yfield = "Pressure";
      auto fun = [&](double rho, double p) {
        double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
        return vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      xfield = "Density";
      yfield = "Temperature";
      auto fun = [&](double rho, double T) {
        double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T);
        return vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::DP_DRHO) {
    zfield = "PressureDerivativeDensity";

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      xfield = "Density";
      yfield = "Specific Internal Energy";
      auto fun = [&](double rho, double e) {return vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      xfield = "Density";
      yfield = "Pressure";
      auto fun = [&](double rho, double p) {
        double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
        return vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      xfield = "Density";
      yfield = "Temperature";
      auto fun = [&](double rho, double T) {
        double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T);
        return vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::BULK_MODULUS) {
    zfield = "Bulk Modulus";

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      xfield = "Density";
      yfield = "Specific Internal Energy";
      auto fun = [&](double rho, double e) {return rho*vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      xfield = "Density";
      yfield = "Pressure";
      auto fun = [&](double rho, double p) {
        double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
        return rho*vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      xfield = "Density";
      yfield = "Temperature";
      auto fun = [&](double rho, double T) {
        double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T);
        return rho*vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::TEMPERATURE) {
    zfield = "Temperature";

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      xfield = "Density";
      yfield = "Specific Internal Energy";
      auto fun = [&](double rho, double e) {return vf[id]->GetTemperature(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      xfield = "Density";
      yfield = "Pressure";
      auto fun = [&](double rho, double p) {
        double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
        return vf[id]->GetTemperature(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::PRESSURE && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      xfield = "Pressure";
      yfield = "Specific Internal Energy";
      auto fun = [&](double p, double e) {
        double rho = vf[id]->GetDensity(p,e);
        return vf[id]->GetTemperature(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }


  if(xfield.empty() || yfield.empty() || zfield.empty()) {
    print_error("*** Error: Unable to create EOS table %d (mat.id: %d). (Try switching the two independent variables.)\n",
                table_id, id);
    exit(-1);
  }

  assert(result.size()>0);

  print("- Outputing EOS Tabulation to %s.\n", tab.filename);
  print("  o Material id: %d, %s as a function of %s and %s.\n", id, zfield.c_str(), xfield.c_str(), yfield.c_str());



  // Open file
  FILE *file = fopen(tab.filename, "w");
  if(!file) {
    print_error("*** Error: Cannot write file %s.\n", tab.filename);
    exit(-1);
  }


  // Extract the file name (without path)
  string filename = tab.filename;
  for(auto loc=filename.end(); loc != filename.begin(); loc--) {
    if(*loc == '/') {
      if(loc+1 != filename.end())
        filename = filename.substr((int)(loc-filename.begin())+1);
      break;
    }
  }



  // Write data to file

  if(result.size()==1) { //1D 
    double vmin, dv;
    if(x0==xmax || Nx==1) {
      print(file, "# Tabulating %s as a function of ", zfield.c_str());
      print(file, "%s [%e] and %s [%e, %e]\n", xfield.c_str(), x0, yfield.c_str(), y0, ymax);
      print(file,"# %s fixed at %e.\n", xfield.c_str(), x0);
      print(file,"# Format: [%s] [%s]\n", yfield.c_str(), zfield.c_str());
      print(file,"# Example Gnuplot commands:\n");
      print(file,"#   plot '%s' using 1:2 with l lw 3 title '%s (%s = %e)'\n", 
            filename.c_str(), zfield.c_str(), xfield.c_str(), x0);
      print(file,"#   set xrange[%e:%e];  set xlabel \"%s\";  set ylabel \"%s\";  set grid;  replot;  replot\n",
            y0, ymax, yfield.c_str(), zfield.c_str());

      vmin = y0;
      dv = Ny==1 ? 0.0 : (ymax-y0)/(Ny-1); 
    } 
    else {
      assert(y0==ymax || Ny==1);
      print(file, "# Tabulating %s as a function of ", zfield.c_str());
      print(file, "%s [%e, %e] and %s [%e]\n", xfield.c_str(), x0, xmax, yfield.c_str(), y0);
      print(file,"# %s fixed at %e.\n", yfield.c_str(), y0);
      print(file,"# Format: [%s] [%s]\n", xfield.c_str(), zfield.c_str());
      print(file,"# Example Gnuplot commands:\n");
      print(file,"#   plot '%s' using 1:2 with l lw 3 title '%s (%s = %e)'\n", 
            filename.c_str(), zfield.c_str(), yfield.c_str(), y0);
      print(file,"#   set xrange[%e:%e];  set xlabel \"%s\";  set ylabel \"%s\";  set grid;  replot;  replot\n",
            x0, xmax, xfield.c_str(), zfield.c_str());

      vmin = x0;
      dv = Nx==1 ? 0.0 : (xmax-x0)/(Nx-1); 
    }
    for(int i=0; i<(int)result[0].size(); i++)
      print(file, "%16.8e  %16.8e\n", vmin+i*dv, result[0][i]);

    fclose(file);

    return;
  }

  
  //2D
  assert(result.size()>1);
  print(file, "# Tabulating %s as a function of ", zfield.c_str());
  print(file, "%s [%e, %e] and %s [%e, %e]\n", xfield.c_str(), x0, xmax, yfield.c_str(), y0, ymax);
  print(file,"# Format: Gnuplot \"Unstructured\" data\n");
  print(file,"# Example Gnuplot commands:\n");
  print(file,"#   splot '%s' matrix nonuniform with pm3d title '%s'\n", filename.c_str(), zfield.c_str());
  print(file,"#   set view map;  set xrange[%e:%e];  set yrange[%e:%e]\n", x0, xmax, y0, ymax);
  print(file,"#   set xlabel \"%s\";  set ylabel \"%s\";  set title \"%s\";  replot\n",
        xfield.c_str(), yfield.c_str(), zfield.c_str()); 
  double dx = Nx<=1 ? 0.0 : (xmax-x0)/(Nx-1);
  double dy = Ny<=1 ? 0.0 : (ymax-y0)/(Ny-1);

  print(file,"%16d", Nx);
  for(int i=0; i<Nx; i++)
    print(file,"  %16e", x0+i*dx); 
  print(file,"\n");

  for(int j=0; j<(int)result.size(); j++) {
    print(file,"%16e", y0+j*dy);
    for(int i=0; i<(int)result[j].size(); i++)
      print(file,"  %16e", result[j][i]); 
    print(file,"\n");
  }

  fclose(file);

}

//------------------------------------------------------------------------------

