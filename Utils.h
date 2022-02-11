#pragma once
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <mpi.h>

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
//! MPI Rank i will print to stdout
//void print(int i, const char format[],...);
void print(MPI_Comm& comm, int i, const char format[],...);
//--------------------------------------------------
//! MPI Rank 0 will print to a file
void print(FILE* fd, const char format[],...);
void print(MPI_Comm& comm, FILE* fd, const char format[],...);
//--------------------------------------------------
//! Check for NAN
template <class T>
inline int m2c_isnan(const T& t) {return (t != t);}
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
//! Give a constant value to an array of known dimension
inline void setValue(double* out, double in, int dim)
{
  for(int i=0; i<dim; i++)
      out[i] = in;
}
//--------------------------------------------------
