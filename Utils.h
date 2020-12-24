#pragma once

#include <stdarg.h>
#include <stdio.h>
#include <string>

using std::string;

/**************************
 * Utility functions
 **************************
*/
//--------------------------------------------------
//! MPI Rank 0 will print to stdout
void print(const char format[],...);
//--------------------------------------------------
//! MPI Rank i will print to stdout
void print(int i, const char format[],...);
//--------------------------------------------------
//! MPI Rank 0 will print to a file
void print(FILE* fd, const char format[],...);
//--------------------------------------------------
//! Check for NAN
template <class T>
inline int m2c_isnan(const T& t) {return (t != t);}
//--------------------------------------------------
//! Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string getCurrentDateTime();
//--------------------------------------------------
//! Print logo and code version
void printLogo();
//--------------------------------------------------
//! Call MPI_Finalize and exit (with error)
void terminate();
