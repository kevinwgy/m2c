/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<PrescribedMotionOperator.h>
#include<SpaceVariable.h>
#include<Vector5D.h>
#include<cassert>
#include<fstream>
#include<sstream>
#include<string>

//----------------------------------------------------------------

PrescribedMotionOperator::PrescribedMotionOperator(ObjectMap<PrescribedMotionData> &pm_data,
                                                   int id_size)
{
  velo.clear();

  print("\n");

  for(auto&& pm : pm_data.dataMap) {

    int myid = pm.second->materialid;

    if(myid<0 || myid>=id_size-1) {//id_size-1 is "INACTIVE_MATERIAL_ID"
      print_error("*** Error: Detected invalid material id (%d) in PrescribedMotion[%d].\n",
                  myid, pm.first);
      exit_mpi();
    }

    if(velo.find(myid) != velo.end()) {
      print_error("*** Error: Detected repeated material id (%d) in PrescribedMotion[%d].\n",
                  myid, pm.first);
      exit_mpi();
    }

    print(" - Setting up prescribed motion for material id %d.\n", myid);

    velo[myid] = VelocityTimeHistory();

    if(strcmp(pm.second->velocity_time_history, "") == 0) { // file not specified
      velo[myid].time.push_back(0.0);

      Vec3D v(pm.second->velocity_x, pm.second->velocity_y, pm.second->velocity_z);
      velo[myid].velocity.push_back(v);
      print("   o Velocity: (%e %e %e)\n", v[0], v[1], v[2]);
    }
    else
      ReadVelocityFromFile(velo[myid], pm.second->velocity_time_history);

  }

}

//----------------------------------------------------------------

PrescribedMotionOperator::~PrescribedMotionOperator()
{ }

//----------------------------------------------------------------

void
PrescribedMotionOperator::Destroy()
{ }

//----------------------------------------------------------------

void
PrescribedMotionOperator::ReadVelocityFromFile(VelocityTimeHistory &vth, const char *filename)
{

  assert(vth.time.empty() && vth.velocity.empty());

  std::ifstream input(filename, std::ifstream::in);
  if(input.fail()) {
    print_error("*** Error: Cannot open velocity time-history file (%s).\n", filename);
    exit_mpi();
  }  

  string line;
  double time;
  Vec3D v;

  int line_number = 0;
  while(getline(input, line)) {
  
    line_number++;

    auto first_nonspace_id = line.find_first_not_of(" ");
    if((unsigned)first_nonspace_id<line.size() && line[first_nonspace_id] == '#')
      continue;  //this line is comment

    std::stringstream linestream(line);

    // read time
    if(!(linestream >> time)) //the first word on the line
      continue;

    // make sure time is in ascending order
    if(vth.time.size()>0 && vth.time.back() >= time) {
      print_error("*** Error: The time points in %s are not in ascending order (%e).\n", filename,
                  time);
      exit_mpi();
    }

    // read velocity vector
    for(int i=0; i<3; i++) {
      if(!(linestream >> v[i])) {
        print_error("*** Error: File %s is broken around line %d.\n", filename, line_number);
        exit_mpi();
      }
    }

    vth.time.push_back(time);
    vth.velocity.push_back(v); 
  }

  if(vth.time.empty()) {
    print_error("*** Error: Did not find any data in file %s.\n", filename);
    exit_mpi();
  }

  print("   o Read velocity time-history from %s (%d time points).\n", filename, vth.time.size());

}

//----------------------------------------------------------------

void
PrescribedMotionOperator::UpdateVelocity(SpaceVariable3D &V, SpaceVariable3D &ID, double time)
{

  int i0, j0, k0, imax, jmax, kmax;
  V.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  Vec5D***  v  = (Vec5D***)V.GetDataPointer();
  double*** id = ID.GetDataPointer();

  int ta, tb;
  double cta, ctb;

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        auto it = velo.find(id[k][j][i]);
        if(it == velo.end())
          continue; // velocity not specified for this id.

        if(it->second.velocity.size() == 1) {//constant velocity 
          for(int p=0; p<3; p++)
            v[k][j][i][1+p] = it->second.velocity[0][p];
        }
        else { // interpolation

          std::vector<double> &ts(it->second.time);
          std::vector<Vec3D> &vs(it->second.velocity);
            
          if(time<=ts.front()) {
            for(int p=0; p<3; p++)
              v[k][j][i][1+p] = vs[0][p];
          }
          else if(time>=ts.back()) {
            for(int p=0; p<3; p++)
              v[k][j][i][1+p] = vs.back()[p];
          }
          else {
            tb = lower_bound(ts.begin(), ts.end(), time) - ts.begin();
            ta = tb - 1;
            cta = (ts[tb]-time)/(ts[tb]-ts[ta]); 
            ctb = 1.0 - cta;
            for(int p=0; p<3; p++)
              v[k][j][i][1+p] = cta*vs[ta][p] + ctb*vs[tb][p];
          }
        }

      }

  ID.RestoreDataPointerToLocalVector(); //no changes made
  V.RestoreDataPointerAndInsert();

}

//----------------------------------------------------------------









//----------------------------------------------------------------


