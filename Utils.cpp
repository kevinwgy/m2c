#include <iostream>
#include <Utils.h>
#include <time.h>
#include <version.h>
#include <stdio.h>
#include <cstring>
#include <cmath> //floor
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
  MPI_Barrier(MPI_COMM_WORLD);
  return;
}

//--------------------------------------------------
// MPI Rank 0 will print to stdout
void print(MPI_Comm& comm, const char format[],...)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(!rank) {
    va_list Argp;
    va_start(Argp, format);
    vprintf(format, Argp);
    va_end(Argp);
  }
  MPI_Barrier(comm);
  return;
}

//--------------------------------------------------
// MPI Rank 0 will print to stdout in red color
void print_error(const char format[],...)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(!rank) {

    char format_colored[strlen(format)+40] = "";
    strcat(format_colored, "\033[0;31m");
    strcat(format_colored, format);
    strcat(format_colored, "\033[0m");

    va_list Argp;
    va_start(Argp, format);
    vprintf(format_colored, Argp);
    va_end(Argp);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return;
}

//--------------------------------------------------
// MPI Rank 0 will print to stdout in red color
void print_error(MPI_Comm& comm, const char format[],...)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  if(!rank) {

    char format_colored[strlen(format)+40] = "";
    strcat(format_colored, "\033[0;31m");
    strcat(format_colored, format);
    strcat(format_colored, "\033[0m");

    va_list Argp;
    va_start(Argp, format);
    vprintf(format_colored, Argp);
    va_end(Argp);
  }
  MPI_Barrier(comm);
  return;
}

//--------------------------------------------------
// MPI Rank i will print to stdout
/*
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

  MPI_Barrier(MPI_COMM_WORLD);
  return;
}
*/
//--------------------------------------------------
// MPI Rank i will print to stdout
void print(MPI_Comm& comm, int i, const char format[],...)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  if(rank == i) {
    va_list Argp;
    va_start(Argp, format);
    vprintf(format, Argp);
    va_end(Argp);
  }

  MPI_Barrier(comm);
  return;
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

  MPI_Barrier(MPI_COMM_WORLD);
  return;
}

//--------------------------------------------------
// MPI Rank 0 will print to a file
void print(MPI_Comm& comm, FILE* fd, const char format[],...)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  if(!rank) {
    va_list Argp;
    va_start(Argp, format);
    vfprintf(fd, format, Argp);
    va_end(Argp);
  }

  MPI_Barrier(comm);
  return;
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
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X ", &tstruct);
    if(strlen(tzname[1]) != 0)
      strcat(buf, tzname[1]); //daylight saving time
    else
      strcat(buf, tzname[0]); //standard time


    return buf;
}

//--------------------------------------------------

//--------------------------------------------------
// Print logo 
void printHeader(int argc, char *argv[])
{
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(!rank) {
    cout << endl;
    cout << "\033[0;36m                                      _..._      \033[0m" << endl;
    cout << "\033[0;36m                      .-''-.       .-'_..._''.   \033[0m" << endl;
    cout << "\033[0;36m   __  __   ___     .' .-.  )    .' .'      '.\\  \033[0m" << endl;
    cout << "\033[0;36m  |  |/  `.'   `.  / .'  / /    / .'             \033[0m" << endl;
    cout << "\033[0;36m  |   .-.  .-.   '(_/   / /    . '               \033[0m" << endl;
    cout << "\033[0;36m  |  |  |  |  |  |     / /     | |               \033[0m" << endl;
    cout << "\033[0;36m  |  |  |  |  |  |    / /      | |               \033[0m" << endl;
    cout << "\033[0;36m  |  |  |  |  |  |   . '       . '               \033[0m" << endl;
    cout << "\033[0;36m  |  |  |  |  |  |  / /    _.-')\\ '.          .  \033[0m" << endl;
    cout << "\033[0;36m  |__|  |__|  |__|.' '  _.'.-''  '. `._____.-'/  \033[0m" << endl;
    cout << "\033[0;36m                 /  /.-'_.'        `-.______ /   \033[0m" << endl;
    cout << "\033[0;36m                /    _.'                    `    \033[0m" << endl;
    cout << "\033[0;36m               ( _.-'                            \033[0m" << endl;
    cout << endl;
    cout << "Revision: " << GIT_REV << " | " << "Branch: " << GIT_BRANCH << " | " << "Tag: " << GIT_TAG << endl;
    cout << "Simulation started at: " << getCurrentDateTime() << endl;
    cout << "Using " << size << " processor cores." << endl;
    cout << "Command:";
    for(int i=0; i<argc; i++)
      cout << " " << argv[i];
    cout << endl;
    cout << endl;
    cout.flush();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return;
}

//--------------------------------------------------
// Terminate program properly
void exit_mpi()
{
  MPI_Finalize();
  exit(-1);
}

//--------------------------------------------------

bool isTimeToWrite(double time, double dt, int time_step, double frequency_dt, int frequency,
                   double last_snapshot_time, bool force_write)
{
  //! If "force_write", just avoid duplication
  if(force_write)
     return (time > last_snapshot_time + 0.1*dt);

  //! Check frequency_dt. If it is not specified, use frequency
  if(frequency_dt > 0) {
    if(floor(time/frequency_dt) != floor((time-dt)/frequency_dt))
      return true;
    else
      return false;
  }
  else { //!< use frequency
    if((frequency > 0) && (time_step % frequency == 0))
      return true;
    else
      return false;
  }
}





