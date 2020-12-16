#include <iostream>
#include <Utils.h>
#include <time.h>
#include <mpi.h>
#include <version.h>
#include <stdio.h>
using std::cout;
using std::endl;
//--------------------------------------------------
// MPI Rank 0 will print to stdout
void print(const char format[],...)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(!rank) {
    va_list Argp;
    va_start(Argp, format);
    vprintf(format, Argp);
    va_end(Argp);
  }
}

//--------------------------------------------------
// MPI Rank i will print to stdout
void print(int i, const char format[],...)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(rank == i) {
    va_list Argp;
    va_start(Argp, format);
    vprintf(format, Argp);
    va_end(Argp);
  }
}

//--------------------------------------------------
// MPI Rank 0 will print to a file
void print(FILE* fd, const char format[],...)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(!rank) {
    va_list Argp;
    va_start(Argp, format);
    vfprintf(fd, format, Argp);
    va_end(Argp);
  }
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

//--------------------------------------------------
// Print logo 
void printLogo()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(!rank) {
    cout << endl;
    cout << "                                     _..._     " << endl;
    cout << "                     .-''-.       .-'_..._''.  " << endl;
    cout << "  __  __   ___     .' .-.  )    .' .'      '.\\  " << endl;
    cout << " |  |/  `.'   `.  / .'  / /    / .'            " << endl;
    cout << " |   .-.  .-.   '(_/   / /    . '              " << endl;
    cout << " |  |  |  |  |  |     / /     | |              " << endl;
    cout << " |  |  |  |  |  |    / /      | |              " << endl;
    cout << " |  |  |  |  |  |   . '       . '              " << endl;
    cout << " |  |  |  |  |  |  / /    _.-')\\ '.          . " << endl;
    cout << " |__|  |__|  |__|.' '  _.'.-''  '. `._____.-'/ " << endl;
    cout << "                /  /.-'_.'        `-.______ /  " << endl;
    cout << "               /    _.'                    `   " << endl;
    cout << "              ( _.-'                           " << endl;
    cout << endl;
    cout << "Revision: " << GIT_REV << " | " << "Branch: " << GIT_BRANCH << " | " << "Tag: " << GIT_TAG << endl;
    cout << "Simulation started at: " << getCurrentDateTime() << endl;
    cout << endl;
    cout.flush();
  }
}

//--------------------------------------------------
