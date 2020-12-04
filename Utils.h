#ifndef _UTILS_H_
#define _UTILS_H_

#include <time.h>
#include <stdarg.h>
#include <stdio.h>
#include <petscsys.h>

//--------------------------------------------------
// MPI Rank 0 will print to stdout
void print(const char format[],...)
{
  va_list Argp;
  va_start(Argp, format);
  PetscPrintf(PETSC_COMM_WORLD, format, Argp);
}

//--------------------------------------------------
// MPI Rank i will print to stdout
void print(int i, const char format[],...)
{
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  va_list Argp;
  va_start(Argp, format);

  if(rank==i)
    PetscPrintf(PETSC_COMM_SELF, format, Argp);
}

//--------------------------------------------------
// All MPI processes will print to the screen
void Print(const char format[],...)
{
  va_list Argp;
  va_start(Argp, format);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, format, Argp); 
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}

//--------------------------------------------------
// MPI Rank 0 will print to a file
void print(FILE* fd, const char format[],...)
{
  va_list Argp;
  va_start(Argp, format);
  PetscFPrintf(PETSC_COMM_WORLD, fd, format, Argp);
}

//--------------------------------------------------
// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string getCurrentDateTime()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}
#endif

