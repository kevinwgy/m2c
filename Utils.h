/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <mpi.h>
#include <cctype>
#include <cassert>

using std::string;

//--------------------------------------------------
/**************************
 * Utility functions
 **************************
*/
//--------------------------------------------------
//! MPI Rank 0 will print to stdout
void print(const char format[],...);
void print(MPI_Comm& comm, const char format[],...);
//--------------------------------------------------
//! MPI Rank 0 will print to stdout in red color
void print_error(const char format[],...);
void print_error(MPI_Comm& comm, const char format[],...);
//--------------------------------------------------
//! MPI Rank 0 will print to stdout in purple color
void print_warning(const char format[],...);
void print_warning(MPI_Comm& comm, const char format[],...);
//--------------------------------------------------
//! MPI Rank i will print to stdout
//void print(int i, const char format[],...);
void print(MPI_Comm& comm, int i, const char format[],...);
//--------------------------------------------------
//! MPI Rank 0 will print to a file
void print(FILE* fd, const char format[],...);
void print(MPI_Comm& comm, FILE* fd, const char format[],...);
//--------------------------------------------------
void mpi_barrier();
//--------------------------------------------------
//! Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string getCurrentDateTime();
//--------------------------------------------------
//! Print logo and code version
void printHeader(int argc, char* argv[]);
//--------------------------------------------------
//! Call MPI_Finalize and exit (with error)
void exit_mpi();
//--------------------------------------------------
//! Determine if solution snapshot should be written
bool isTimeToWrite(double time, double dt, int time_step, double frequency_dt, int frequency,
                   double last_snapshot_time, bool force_write);
//--------------------------------------------------
//! case-insensitive string compare (true: equal;  false: unequal)
inline bool same_strings_insensitive(std::string str1, std::string str2)
{
  return ((str1.size() == str2.size()) && 
           std::equal(str1.begin(), str1.end(), str2.begin(), 
                      [](char &c1, char &c2){return (c1==c2 || std::toupper(c1)==std::toupper(c2));}));
}
//--------------------------------------------------
//! Copy a double array of known dimension
inline void copyarray(double* in, double* out, int dim)
{
  for(int i=0; i<dim; i++)
      out[i] = in[i];
}
//--------------------------------------------------
//! Add a double array of known dimension
inline void addarray(double* in, double* out, int dim)
{
  for(int i=0; i<dim; i++)
      out[i] += in[i];
}
//--------------------------------------------------
//! Copy a double array of known dimension, flip the sign of one element
inline void copyarray_flip(double* in, double* out, int dim, int flipdim)
{
  for(int i=0; i<dim; i++)
    out[i] = (i == flipdim) ? -in[i] : in[i];
}
//--------------------------------------------------
//! Copy a double array of known dimension, flip the sign of multiple consecutive elements 
inline void copyarray_flip(double* in, double* out, int dim, int flipdim0, int flipsize)
{
  for(int i=0; i<dim; i++) 
    out[i] = (i>=flipdim0 && i<flipdim0+flipsize) ? -in[i] : in[i];
}
//--------------------------------------------------
//! Give a constant value to an array of known dimension
inline void setValue(double* out, double in, int dim)
{
  for(int i=0; i<dim; i++)
      out[i] = in;
}
//--------------------------------------------------

template<typename Functor>
void tabulate2Dfunction_uniform(Functor fun, double xmin, double xmax,
                                int Nx, double ymin, double ymax, int Ny,
                                std::vector<std::vector<double> > &result)
{
  //NOTE: fun must be of type "double fun(double x, double y)".

  //Check for the case of a 1D plot
  if(xmin==xmax || Nx==1) {
    result.resize(1);
    assert(Ny>=1);
    double dy = Ny==1 ? 0.0 : (ymax-ymin)/(Ny-1);
    result[0].reserve(Ny);
    for(int i=0; i<Ny; i++)
      result[0].push_back(fun(xmin, ymin+i*dy));
    return;
  }
  else if(ymin==ymax || Ny==1) {
    result.resize(1);
    assert(Nx>=1);
    double dx = Nx==1 ? 0.0 : (xmax-xmin)/(Nx-1);
    result[0].reserve(Nx);
    for(int i=0; i<Nx; i++)
      result[0].push_back(fun(xmin+i*dx, ymin));
    return;
  }

  // 2D
  assert(Nx>1 && Ny>1);
  double dx = (xmax-xmin)/(Nx-1);
  double dy = (ymax-ymin)/(Ny-1);
  result.resize(Ny);
  double y = ymin;
  for(auto&& res : result) {
    res.reserve(Nx);
    for(int i=0; i<Nx; i++)
      res.push_back(fun(xmin+i*dx, y));
    y += dy;
  }
}

//--------------------------------------------------


//--------------------------------------------------





