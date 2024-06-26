/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<TerminalVisualization.h>
#include<string.h>
#include<Vector5D.h>
#include<bits/stdc++.h>
#include<unistd.h> //sleep

using std::vector;

//------------------------------------------------------------------

TerminalVisualization::TerminalVisualization(MPI_Comm &comm_, TerminalVisualizationData &iod_terminal_, 
                                             GlobalMeshInfo &global_mesh_, vector<VarFcnBase*> &vf_, 
                                             IonizationOperator* ion_)
                     : comm(comm_), iod_terminal(iod_terminal_), global_mesh(global_mesh_), vf(vf_), ion(ion_),
                       myout(NULL)
{
  
  if(iod_terminal.plane == TerminalVisualizationData::NONE)
    return; //not activated

  if(strcmp(iod_terminal.filename, "")) {
    print(comm, "\n- Initializing terminal visualization in %s.\n", iod_terminal.filename);
    print(comm, "  o Command for visualization: 'tail -f %s'\n"
                "    or 'less -R %s' (then, press shift+f)\n", iod_terminal.filename, iod_terminal.filename);
  } else
    print(comm, "\n- Initializing terminal visualization.\n");


  // check for user errors
  if(iod_terminal.pause<0.0) {
    print_error(comm, "*** Error: 'Pause' can not be set to a negative number (%e).\n", iod_terminal.pause);
    exit(-1);
  }
  if(iod_terminal.variable == TerminalVisualizationData::MEANCHARGE && !ion) {
    print_error(comm, "*** Error: Cannot print mean charge number. Ionization solver is not activated.\n");
    exit(-1);
  }

  // get mpi info
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  iFrame = 0;
  last_snapshot_time = -1.0;
  last_snapshot_clocktime = -1.0;

  // determine whether use the clock time or the simulation/physical time
  if(iod_terminal.frequency<=0.0 && iod_terminal.frequency_dt<=0.0) {
    if(iod_terminal.frequency_clocktime<=0) {
      print_error(comm, "*** Error: Frequency of terminal visualization is not specified properly.\n");
      exit(-1);
    }
    use_clocktime = true;
  } else
    use_clocktime = false; 

  // determine the region to visualize
  double coord(0.0);
  double hmin(0.0), hmax(0.0), vmin(0.0), vmax(0.0);
  int hindex(-1), vindex(-1), pindex(-1); //0,1,2
  switch (iod_terminal.plane) {
    case TerminalVisualizationData::YZ :
      if(iod_terminal.coordinate < global_mesh.x_glob.front() || 
         iod_terminal.coordinate > global_mesh.x_glob.back()) 
        coord = 0.5*(global_mesh.x_glob.front() + global_mesh.x_glob.back());
      else
        coord = iod_terminal.coordinate;

      hmin = std::max(iod_terminal.horizontal_min, global_mesh.y_glob.front());
      hmax = std::min(iod_terminal.horizontal_max, global_mesh.y_glob.back());
      vmin = std::max(iod_terminal.vertical_min,   global_mesh.z_glob.front());
      vmax = std::min(iod_terminal.vertical_max,   global_mesh.z_glob.back());

      // these two are member variables of the class
      xyzmin = Vec3D(coord, hmin, vmin);
      xyzmax = Vec3D(coord, hmax, vmax);

      hindex = 1; //y
      vindex = 2; //z
      pindex = 0; //x

      break;

    case TerminalVisualizationData::XZ :
      if(iod_terminal.coordinate < global_mesh.y_glob.front() || 
         iod_terminal.coordinate > global_mesh.y_glob.back()) 
        coord = 0.5*(global_mesh.y_glob.front() + global_mesh.y_glob.back());
      else
        coord = iod_terminal.coordinate;

      hmin = std::max(iod_terminal.horizontal_min, global_mesh.x_glob.front());
      hmax = std::min(iod_terminal.horizontal_max, global_mesh.x_glob.back());
      vmin = std::max(iod_terminal.vertical_min,   global_mesh.z_glob.front());
      vmax = std::min(iod_terminal.vertical_max,   global_mesh.z_glob.back());

      // these two are member variables of the class
      xyzmin = Vec3D(hmin, coord, vmin);
      xyzmax = Vec3D(hmax, coord, vmax);

      hindex = 0; //x
      vindex = 2; //z
      pindex = 1; //y

      break;

    case TerminalVisualizationData::XY :
      if(iod_terminal.coordinate < global_mesh.z_glob.front() || 
         iod_terminal.coordinate > global_mesh.z_glob.back()) 
        coord = 0.5*(global_mesh.z_glob.front() + global_mesh.z_glob.back());
      else
        coord = iod_terminal.coordinate;

      hmin = std::max(iod_terminal.horizontal_min, global_mesh.x_glob.front());
      hmax = std::min(iod_terminal.horizontal_max, global_mesh.x_glob.back());
      vmin = std::max(iod_terminal.vertical_min,   global_mesh.y_glob.front());
      vmax = std::min(iod_terminal.vertical_max,   global_mesh.y_glob.back());

      // these two are member variables of the class
      xyzmin = Vec3D(hmin, vmin, coord);
      xyzmax = Vec3D(hmax, vmax, coord);

      hindex = 0; //x
      vindex = 1; //y
      pindex = 2; //z

      break;

    case TerminalVisualizationData::NONE : //just to avoid compiler warning
      break;
  }

  //default dimensions 
  int nmax = 128;
  int nmin = 32; 

  //calculate dx (same for y and z);
  dx = iod_terminal.dx;
  for(int i=0; i<3; i++) {
    if(xyzmin[i] == xyzmax[i])
      continue;
    double dx_min = (xyzmax[i] - xyzmin[i])/(nmax-1);
    double dx_max = (xyzmax[i] - xyzmin[i])/(nmin-1);
    if(dx<dx_min)
      dx = dx_min;
    else if(dx>dx_max)
      dx = dx_max;
  }
  if(iod_terminal.dx>0 && iod_terminal.dx != dx)
    print(comm, "  o Adjusted the resolution for terminal visualization. dx: %e -> %e.\n", 
          iod_terminal.dx, dx);

  // determine nrows and ncols
  nrows = (int)((vmax - vmin)/dx) + 1;
  ncols = (int)((hmax - hmin)/dx) + 1;
  dv = (vmax - vmin)/(nrows-1);
  dh = (hmax - hmin)/(ncols-1); //dv and dh should be similar, and similar to dx

  // determine ijk
  ijk.resize(nrows*ncols);
  Vec3D xyz;
  xyz[pindex] = coord;
  for(int j=0; j<nrows; j++) {
    xyz[vindex] = vmax - j*dv;
    for(int i=0; i<ncols; i++) {
      xyz[hindex] = hmin + i*dh;
      // find closest node to xyz
      ijk[j*ncols+i] = global_mesh.FindClosestNodeToPoint(xyz, false); //not including ghosts
    }
  }


  // setup colormap
  if(iod_terminal.colormap == TerminalVisualizationData::GRAYSCALE)
    SetupGrayColorMap();
  else
    SetupTurboColorMap();



  // create ostream
  std::streambuf *buf;
  if(strcmp(iod_terminal.filename, "")) {
    outfile.open(iod_terminal.filename, std::ofstream::out);
    if(!outfile.is_open()) {
      fprintf(stdout,"\033[0;31m*** Error: Unable to open file %s for visualization.\n", iod_terminal.filename);
      exit(-1);
    }
    buf = outfile.rdbuf();
  } else
    buf = std::cout.rdbuf();

  myout = new std::ostream(buf);

}

//------------------------------------------------------------------

TerminalVisualization::~TerminalVisualization()
{
  if(outfile.is_open())
    outfile.close();  
  if(myout) {
    delete myout;
  }
}

//------------------------------------------------------------------
// The 5; (256 color) format starts with the 16 original colors as colors 0-15. 
// Color numbers 0 to 7 are the default terminal colors, the actual RGB value of 
// which is not standardized and can often be configured. Color numbers 8 to 15 
// are the "bright" colors. Most of the time these are a lighter shade of the 
// color with index - 8. They are also not standardized and can often be configured. 
// Depending on terminal and shell, they are often used instead of or in conjunction 
// with bold font faces.
// Color numbers 16 to 231 are RGB colors. These 216 colors are defined by 6 values 
// on each of the three RGB axes. That is, instead of values 0 - 255, each color only 
// ranges from 0 - 5. The color number is then calculated like this:
// number = 16 + 36 * r + 6 * g + b, with r, g and b in the range 0 - 5.
// The color numbers 232 to 255 are grayscale with 24 shades of gray from dark to light.
//------------------------------------------------------------------

void
TerminalVisualization::SetupGrayColorMap()
{
  number_colors = 24;
  gray_ANSI.reserve(number_colors);

  for(int i=0; i<number_colors; i++)
    gray_ANSI.push_back(232+i);

  under_color_ANSI = 0; //black
  over_color_ANSI  = 15; //white
}

//------------------------------------------------------------------

void
TerminalVisualization::SetupTurboColorMap()
{
  SetupTurboRGB(); //build the turbo_rbg vector

  number_colors = 24;
  turbo_ANSI.reserve(number_colors);

  double step = turbo_rgb.size()/number_colors; //this is "double"
  int r,g,b;
  for(int i=0; i<number_colors; i++) {
    Vec3D& rgb(turbo_rgb[round((i+0.5)*step)]);
    // map each color (r,g,b) from (0-1) -> (0,1,2,3,4,5)
    r = (int)(rgb[0]*6 - 1e-8); //to avoid r = 6
    g = (int)(rgb[1]*6 - 1e-8);
    b = (int)(rgb[2]*6 - 1e-8);
    turbo_ANSI.push_back(16 + 36*r + 6*g + b);
  }

  under_color_ANSI = 0; //black
  over_color_ANSI  = 15; //white
}

//------------------------------------------------------------------

void
TerminalVisualization::SetupTurboRGB()
{

  // Got from internet. We actually don't need this many points here.
  turbo_rgb = {
    Vec3D(0.18995, 0.07176, 0.23217),
    Vec3D(0.19483, 0.08339, 0.26149),
    Vec3D(0.19956, 0.09498, 0.29024),
    Vec3D(0.20415, 0.10652, 0.31844),
    Vec3D(0.20860, 0.11802, 0.34607),
    Vec3D(0.21291, 0.12947, 0.37314),
    Vec3D(0.21708, 0.14087, 0.39964),
    Vec3D(0.22111, 0.15223, 0.42558),
    Vec3D(0.22500, 0.16354, 0.45096),
    Vec3D(0.22875, 0.17481, 0.47578),
    Vec3D(0.23236, 0.18603, 0.50004),
    Vec3D(0.23582, 0.19720, 0.52373),
    Vec3D(0.23915, 0.20833, 0.54686),
    Vec3D(0.24234, 0.21941, 0.56942),
    Vec3D(0.24539, 0.23044, 0.59142),
    Vec3D(0.24830, 0.24143, 0.61286),
    Vec3D(0.25107, 0.25237, 0.63374),
    Vec3D(0.25369, 0.26327, 0.65406),
    Vec3D(0.25618, 0.27412, 0.67381),
    Vec3D(0.25853, 0.28492, 0.69300),
    Vec3D(0.26074, 0.29568, 0.71162),
    Vec3D(0.26280, 0.30639, 0.72968),
    Vec3D(0.26473, 0.31706, 0.74718),
    Vec3D(0.26652, 0.32768, 0.76412),
    Vec3D(0.26816, 0.33825, 0.78050),
    Vec3D(0.26967, 0.34878, 0.79631),
    Vec3D(0.27103, 0.35926, 0.81156),
    Vec3D(0.27226, 0.36970, 0.82624),
    Vec3D(0.27334, 0.38008, 0.84037),
    Vec3D(0.27429, 0.39043, 0.85393),
    Vec3D(0.27509, 0.40072, 0.86692),
    Vec3D(0.27576, 0.41097, 0.87936),
    Vec3D(0.27628, 0.42118, 0.89123),
    Vec3D(0.27667, 0.43134, 0.90254),
    Vec3D(0.27691, 0.44145, 0.91328),
    Vec3D(0.27701, 0.45152, 0.92347),
    Vec3D(0.27698, 0.46153, 0.93309),
    Vec3D(0.27680, 0.47151, 0.94214),
    Vec3D(0.27648, 0.48144, 0.95064),
    Vec3D(0.27603, 0.49132, 0.95857),
    Vec3D(0.27543, 0.50115, 0.96594),
    Vec3D(0.27469, 0.51094, 0.97275),
    Vec3D(0.27381, 0.52069, 0.97899),
    Vec3D(0.27273, 0.53040, 0.98461),
    Vec3D(0.27106, 0.54015, 0.98930),
    Vec3D(0.26878, 0.54995, 0.99303),
    Vec3D(0.26592, 0.55979, 0.99583),
    Vec3D(0.26252, 0.56967, 0.99773),
    Vec3D(0.25862, 0.57958, 0.99876),
    Vec3D(0.25425, 0.58950, 0.99896),
    Vec3D(0.24946, 0.59943, 0.99835),
    Vec3D(0.24427, 0.60937, 0.99697),
    Vec3D(0.23874, 0.61931, 0.99485),
    Vec3D(0.23288, 0.62923, 0.99202),
    Vec3D(0.22676, 0.63913, 0.98851),
    Vec3D(0.22039, 0.64901, 0.98436),
    Vec3D(0.21382, 0.65886, 0.97959),
    Vec3D(0.20708, 0.66866, 0.97423),
    Vec3D(0.20021, 0.67842, 0.96833),
    Vec3D(0.19326, 0.68812, 0.96190),
    Vec3D(0.18625, 0.69775, 0.95498),
    Vec3D(0.17923, 0.70732, 0.94761),
    Vec3D(0.17223, 0.71680, 0.93981),
    Vec3D(0.16529, 0.72620, 0.93161),
    Vec3D(0.15844, 0.73551, 0.92305),
    Vec3D(0.15173, 0.74472, 0.91416),
    Vec3D(0.14519, 0.75381, 0.90496),
    Vec3D(0.13886, 0.76279, 0.89550),
    Vec3D(0.13278, 0.77165, 0.88580),
    Vec3D(0.12698, 0.78037, 0.87590),
    Vec3D(0.12151, 0.78896, 0.86581),
    Vec3D(0.11639, 0.79740, 0.85559),
    Vec3D(0.11167, 0.80569, 0.84525),
    Vec3D(0.10738, 0.81381, 0.83484),
    Vec3D(0.10357, 0.82177, 0.82437),
    Vec3D(0.10026, 0.82955, 0.81389),
    Vec3D(0.09750, 0.83714, 0.80342),
    Vec3D(0.09532, 0.84455, 0.79299),
    Vec3D(0.09377, 0.85175, 0.78264),
    Vec3D(0.09287, 0.85875, 0.77240),
    Vec3D(0.09267, 0.86554, 0.76230),
    Vec3D(0.09320, 0.87211, 0.75237),
    Vec3D(0.09451, 0.87844, 0.74265),
    Vec3D(0.09662, 0.88454, 0.73316),
    Vec3D(0.09958, 0.89040, 0.72393),
    Vec3D(0.10342, 0.89600, 0.71500),
    Vec3D(0.10815, 0.90142, 0.70599),
    Vec3D(0.11374, 0.90673, 0.69651),
    Vec3D(0.12014, 0.91193, 0.68660),
    Vec3D(0.12733, 0.91701, 0.67627),
    Vec3D(0.13526, 0.92197, 0.66556),
    Vec3D(0.14391, 0.92680, 0.65448),
    Vec3D(0.15323, 0.93151, 0.64308),
    Vec3D(0.16319, 0.93609, 0.63137),
    Vec3D(0.17377, 0.94053, 0.61938),
    Vec3D(0.18491, 0.94484, 0.60713),
    Vec3D(0.19659, 0.94901, 0.59466),
    Vec3D(0.20877, 0.95304, 0.58199),
    Vec3D(0.22142, 0.95692, 0.56914),
    Vec3D(0.23449, 0.96065, 0.55614),
    Vec3D(0.24797, 0.96423, 0.54303),
    Vec3D(0.26180, 0.96765, 0.52981),
    Vec3D(0.27597, 0.97092, 0.51653),
    Vec3D(0.29042, 0.97403, 0.50321),
    Vec3D(0.30513, 0.97697, 0.48987),
    Vec3D(0.32006, 0.97974, 0.47654),
    Vec3D(0.33517, 0.98234, 0.46325),
    Vec3D(0.35043, 0.98477, 0.45002),
    Vec3D(0.36581, 0.98702, 0.43688),
    Vec3D(0.38127, 0.98909, 0.42386),
    Vec3D(0.39678, 0.99098, 0.41098),
    Vec3D(0.41229, 0.99268, 0.39826),
    Vec3D(0.42778, 0.99419, 0.38575),
    Vec3D(0.44321, 0.99551, 0.37345),
    Vec3D(0.45854, 0.99663, 0.36140),
    Vec3D(0.47375, 0.99755, 0.34963),
    Vec3D(0.48879, 0.99828, 0.33816),
    Vec3D(0.50362, 0.99879, 0.32701),
    Vec3D(0.51822, 0.99910, 0.31622),
    Vec3D(0.53255, 0.99919, 0.30581),
    Vec3D(0.54658, 0.99907, 0.29581),
    Vec3D(0.56026, 0.99873, 0.28623),
    Vec3D(0.57357, 0.99817, 0.27712),
    Vec3D(0.58646, 0.99739, 0.26849),
    Vec3D(0.59891, 0.99638, 0.26038),
    Vec3D(0.61088, 0.99514, 0.25280),
    Vec3D(0.62233, 0.99366, 0.24579),
    Vec3D(0.63323, 0.99195, 0.23937),
    Vec3D(0.64362, 0.98999, 0.23356),
    Vec3D(0.65394, 0.98775, 0.22835),
    Vec3D(0.66428, 0.98524, 0.22370),
    Vec3D(0.67462, 0.98246, 0.21960),
    Vec3D(0.68494, 0.97941, 0.21602),
    Vec3D(0.69525, 0.97610, 0.21294),
    Vec3D(0.70553, 0.97255, 0.21032),
    Vec3D(0.71577, 0.96875, 0.20815),
    Vec3D(0.72596, 0.96470, 0.20640),
    Vec3D(0.73610, 0.96043, 0.20504),
    Vec3D(0.74617, 0.95593, 0.20406),
    Vec3D(0.75617, 0.95121, 0.20343),
    Vec3D(0.76608, 0.94627, 0.20311),
    Vec3D(0.77591, 0.94113, 0.20310),
    Vec3D(0.78563, 0.93579, 0.20336),
    Vec3D(0.79524, 0.93025, 0.20386),
    Vec3D(0.80473, 0.92452, 0.20459),
    Vec3D(0.81410, 0.91861, 0.20552),
    Vec3D(0.82333, 0.91253, 0.20663),
    Vec3D(0.83241, 0.90627, 0.20788),
    Vec3D(0.84133, 0.89986, 0.20926),
    Vec3D(0.85010, 0.89328, 0.21074),
    Vec3D(0.85868, 0.88655, 0.21230),
    Vec3D(0.86709, 0.87968, 0.21391),
    Vec3D(0.87530, 0.87267, 0.21555),
    Vec3D(0.88331, 0.86553, 0.21719),
    Vec3D(0.89112, 0.85826, 0.21880),
    Vec3D(0.89870, 0.85087, 0.22038),
    Vec3D(0.90605, 0.84337, 0.22188),
    Vec3D(0.91317, 0.83576, 0.22328),
    Vec3D(0.92004, 0.82806, 0.22456),
    Vec3D(0.92666, 0.82025, 0.22570),
    Vec3D(0.93301, 0.81236, 0.22667),
    Vec3D(0.93909, 0.80439, 0.22744),
    Vec3D(0.94489, 0.79634, 0.22800),
    Vec3D(0.95039, 0.78823, 0.22831),
    Vec3D(0.95560, 0.78005, 0.22836),
    Vec3D(0.96049, 0.77181, 0.22811),
    Vec3D(0.96507, 0.76352, 0.22754),
    Vec3D(0.96931, 0.75519, 0.22663),
    Vec3D(0.97323, 0.74682, 0.22536),
    Vec3D(0.97679, 0.73842, 0.22369),
    Vec3D(0.98000, 0.73000, 0.22161),
    Vec3D(0.98289, 0.72140, 0.21918),
    Vec3D(0.98549, 0.71250, 0.21650),
    Vec3D(0.98781, 0.70330, 0.21358),
    Vec3D(0.98986, 0.69382, 0.21043),
    Vec3D(0.99163, 0.68408, 0.20706),
    Vec3D(0.99314, 0.67408, 0.20348),
    Vec3D(0.99438, 0.66386, 0.19971),
    Vec3D(0.99535, 0.65341, 0.19577),
    Vec3D(0.99607, 0.64277, 0.19165),
    Vec3D(0.99654, 0.63193, 0.18738),
    Vec3D(0.99675, 0.62093, 0.18297),
    Vec3D(0.99672, 0.60977, 0.17842),
    Vec3D(0.99644, 0.59846, 0.17376),
    Vec3D(0.99593, 0.58703, 0.16899),
    Vec3D(0.99517, 0.57549, 0.16412),
    Vec3D(0.99419, 0.56386, 0.15918),
    Vec3D(0.99297, 0.55214, 0.15417),
    Vec3D(0.99153, 0.54036, 0.14910),
    Vec3D(0.98987, 0.52854, 0.14398),
    Vec3D(0.98799, 0.51667, 0.13883),
    Vec3D(0.98590, 0.50479, 0.13367),
    Vec3D(0.98360, 0.49291, 0.12849),
    Vec3D(0.98108, 0.48104, 0.12332),
    Vec3D(0.97837, 0.46920, 0.11817),
    Vec3D(0.97545, 0.45740, 0.11305),
    Vec3D(0.97234, 0.44565, 0.10797),
    Vec3D(0.96904, 0.43399, 0.10294),
    Vec3D(0.96555, 0.42241, 0.09798),
    Vec3D(0.96187, 0.41093, 0.09310),
    Vec3D(0.95801, 0.39958, 0.08831),
    Vec3D(0.95398, 0.38836, 0.08362),
    Vec3D(0.94977, 0.37729, 0.07905),
    Vec3D(0.94538, 0.36638, 0.07461),
    Vec3D(0.94084, 0.35566, 0.07031),
    Vec3D(0.93612, 0.34513, 0.06616),
    Vec3D(0.93125, 0.33482, 0.06218),
    Vec3D(0.92623, 0.32473, 0.05837),
    Vec3D(0.92105, 0.31489, 0.05475),
    Vec3D(0.91572, 0.30530, 0.05134),
    Vec3D(0.91024, 0.29599, 0.04814),
    Vec3D(0.90463, 0.28696, 0.04516),
    Vec3D(0.89888, 0.27824, 0.04243),
    Vec3D(0.89298, 0.26981, 0.03993),
    Vec3D(0.88691, 0.26152, 0.03753),
    Vec3D(0.88066, 0.25334, 0.03521),
    Vec3D(0.87422, 0.24526, 0.03297),
    Vec3D(0.86760, 0.23730, 0.03082),
    Vec3D(0.86079, 0.22945, 0.02875),
    Vec3D(0.85380, 0.22170, 0.02677),
    Vec3D(0.84662, 0.21407, 0.02487),
    Vec3D(0.83926, 0.20654, 0.02305),
    Vec3D(0.83172, 0.19912, 0.02131),
    Vec3D(0.82399, 0.19182, 0.01966),
    Vec3D(0.81608, 0.18462, 0.01809),
    Vec3D(0.80799, 0.17753, 0.01660),
    Vec3D(0.79971, 0.17055, 0.01520),
    Vec3D(0.79125, 0.16368, 0.01387),
    Vec3D(0.78260, 0.15693, 0.01264),
    Vec3D(0.77377, 0.15028, 0.01148),
    Vec3D(0.76476, 0.14374, 0.01041),
    Vec3D(0.75556, 0.13731, 0.00942),
    Vec3D(0.74617, 0.13098, 0.00851),
    Vec3D(0.73661, 0.12477, 0.00769),
    Vec3D(0.72686, 0.11867, 0.00695),
    Vec3D(0.71692, 0.11268, 0.00629),
    Vec3D(0.70680, 0.10680, 0.00571),
    Vec3D(0.69650, 0.10102, 0.00522),
    Vec3D(0.68602, 0.09536, 0.00481),
    Vec3D(0.67535, 0.08980, 0.00449),
    Vec3D(0.66449, 0.08436, 0.00424),
    Vec3D(0.65345, 0.07902, 0.00408),
    Vec3D(0.64223, 0.07380, 0.00401),
    Vec3D(0.63082, 0.06868, 0.00401),
    Vec3D(0.61923, 0.06367, 0.00410),
    Vec3D(0.60746, 0.05878, 0.00427),
    Vec3D(0.59550, 0.05399, 0.00453),
    Vec3D(0.58336, 0.04931, 0.00486),
    Vec3D(0.57103, 0.04474, 0.00529),
    Vec3D(0.55852, 0.04028, 0.00579),
    Vec3D(0.54583, 0.03593, 0.00638),
    Vec3D(0.53295, 0.03169, 0.00705),
    Vec3D(0.51989, 0.02756, 0.00780),
    Vec3D(0.50664, 0.02354, 0.00863),
    Vec3D(0.49321, 0.01963, 0.00955),
    Vec3D(0.47960, 0.01583, 0.01055)
  };

}

//------------------------------------------------------------------

void
TerminalVisualization::PrintSolutionSnapshot(double time, double dt, int time_step, SpaceVariable3D &V, 
                                             SpaceVariable3D &ID, vector<SpaceVariable3D*> &Phi, 
                                             SpaceVariable3D *L, bool force_write)
{

  if(iod_terminal.plane == TerminalVisualizationData::NONE) //not activated
    return;

  // do we need to print at this time?
  if(!force_write) {
    if(use_clocktime) {
      double current_time = walltime();
      MPI_Allreduce(MPI_IN_PLACE, &current_time, 1, MPI_DOUBLE, MPI_MAX, comm); 
      if(current_time-last_snapshot_clocktime < iod_terminal.frequency_clocktime)
        return;
    } else {
      if(!isTimeToWrite(time, dt, time_step, iod_terminal.frequency_dt, iod_terminal.frequency, 
                        last_snapshot_time, force_write))
        return;
    }
  }


  string solution_name;
  sol.assign(nrows*ncols, 0.0);

  if(iod_terminal.variable == TerminalVisualizationData::DENSITY) {
    int i,j,k;
    solution_name = "Density";
    Vec5D*** v = (Vec5D***)V.GetDataPointer();
    for(int n=0; n<(int)ijk.size(); n++) {
      i = ijk[n][0];
      j = ijk[n][1];
      k = ijk[n][2];
      if(ID.IsHere(i,j,k,false))
        sol[n] = v[k][j][i][0];
    }
    V.RestoreDataPointerToLocalVector();
  }
  else if(iod_terminal.variable == TerminalVisualizationData::VELOCITY) {
    int i,j,k;
    solution_name = "Velocity (Magnitude)";
    Vec5D*** v = (Vec5D***)V.GetDataPointer();
    for(int n=0; n<(int)ijk.size(); n++) {
      i = ijk[n][0];
      j = ijk[n][1];
      k = ijk[n][2];
      if(ID.IsHere(i,j,k,false))
        sol[n] = sqrt(v[k][j][i][1]*v[k][j][i][1] + v[k][j][i][2]*v[k][j][i][2] +
                      v[k][j][i][3]*v[k][j][i][3]);
    }
    V.RestoreDataPointerToLocalVector();
  }
  else if(iod_terminal.variable == TerminalVisualizationData::PRESSURE) {
    int i,j,k;
    solution_name = "Pressure";
    Vec5D*** v = (Vec5D***)V.GetDataPointer();
    for(int n=0; n<(int)ijk.size(); n++) {
      i = ijk[n][0];
      j = ijk[n][1];
      k = ijk[n][2];
      if(ID.IsHere(i,j,k,false))
        sol[n] = v[k][j][i][4]; 
    }
    V.RestoreDataPointerToLocalVector();
  }
  else if(iod_terminal.variable == TerminalVisualizationData::TEMPERATURE) {
    int i,j,k;
    solution_name = "Temperature";
    Vec5D*** v = (Vec5D***)V.GetDataPointer();
    double*** id = ID.GetDataPointer();
    for(int n=0; n<(int)ijk.size(); n++) {
      i = ijk[n][0];
      j = ijk[n][1];
      k = ijk[n][2];
      if(ID.IsHere(i,j,k,false)) {
        double e = vf[id[k][j][i]]->GetInternalEnergyPerUnitMass(v[k][j][i][0],
                                            v[k][j][i][4]);
        sol[n] = vf[id[k][j][i]]->GetTemperature(v[k][j][i][0], e);
      }
    }
    V.RestoreDataPointerToLocalVector();
    ID.RestoreDataPointerToLocalVector();
  }
  else if(iod_terminal.variable == TerminalVisualizationData::MATERIALID) {
    int i,j,k;
    solution_name = "Material ID";
    double*** id = ID.GetDataPointer();
    for(int n=0; n<(int)ijk.size(); n++) {
      i = ijk[n][0];
      j = ijk[n][1];
      k = ijk[n][2];
      if(ID.IsHere(i,j,k,false))
        sol[n] = id[k][j][i];
    }
    ID.RestoreDataPointerToLocalVector();
  }
  else if(iod_terminal.variable == TerminalVisualizationData::LASERRADIANCE) {
    int i,j,k;
    solution_name = "Laser Radiance";
    if(!L) {
      print_error("*** Error: Unable to visualize laser radiance (not activated).\n");
      exit(-1);
    }
    double*** laser = L->GetDataPointer();
    for(int n=0; n<(int)ijk.size(); n++) {
      i = ijk[n][0];
      j = ijk[n][1];
      k = ijk[n][2];
      if(ID.IsHere(i,j,k,false))
        sol[n] = laser[k][j][i];
    }
    L->RestoreDataPointerToLocalVector();
  }
  else if(iod_terminal.variable == TerminalVisualizationData::LEVELSET0) {
    int i,j,k;
    solution_name = "Level Set (0)";
    if(Phi.size()<1) {
      print_error("*** Error: Unable to visualize level set phi[0] (not activated).\n");
      exit(-1);
    }
    double*** phi = Phi[0]->GetDataPointer();
    for(int n=0; n<(int)ijk.size(); n++) {
      i = ijk[n][0];
      j = ijk[n][1];
      k = ijk[n][2];
      if(ID.IsHere(i,j,k,false))
        sol[n] = phi[k][j][i];
    }
    Phi[0]->RestoreDataPointerToLocalVector();
  }
  else if(iod_terminal.variable == TerminalVisualizationData::LEVELSET1) {
    int i,j,k;
    solution_name = "Level Set (1)";
    if(Phi.size()<2) {
      print_error("*** Error: Unable to visualize level set phi[1] (not activated).\n");
      exit(-1);
    }
    double*** phi = Phi[1]->GetDataPointer();
    for(int n=0; n<(int)ijk.size(); n++) {
      i = ijk[n][0];
      j = ijk[n][1];
      k = ijk[n][2];
      if(ID.IsHere(i,j,k,false))
        sol[n] = phi[k][j][i];
    }
    Phi[1]->RestoreDataPointerToLocalVector();
  }
  else if(iod_terminal.variable == TerminalVisualizationData::MEANCHARGE) {
    int i,j,k;
    solution_name = "Mean Charge Number";
    if(!ion) {
      print_error("*** Error: Unable to visualize mean charge number (not activated).\n");
      exit(-1);
    }
    Vec5D*** v = (Vec5D***)V.GetDataPointer();
    double*** id = ID.GetDataPointer();
    for(int n=0; n<(int)ijk.size(); n++) {
      i = ijk[n][0];
      j = ijk[n][1];
      k = ijk[n][2];
      if(ID.IsHere(i,j,k,false)) 
        sol[n] = (ion->ComputeIonizationAtOnePoint((int)id[k][j][i], v[k][j][i][0], v[k][j][i][4]))[0];
    }
    V.RestoreDataPointerToLocalVector();
    ID.RestoreDataPointerToLocalVector();
  }

  if(mpi_rank==0)
    MPI_Reduce(MPI_IN_PLACE, (double*)sol.data(), sol.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
  else
    MPI_Reduce((double*)sol.data(), NULL, sol.size(), MPI_DOUBLE, MPI_SUM, 0, comm);


  // proc #0 prints
  if(!mpi_rank) {

    // determine the scale
    vector<double> tmp = sol;
    std::sort(tmp.begin(), tmp.end());

    double smin = tmp.front();
    double smax = tmp.back();
    double s01  = tmp[2];
    double s99  = tmp[tmp.size()-3];
    double range = s99 - s01;
    bool below_min(false), above_max(false);
    if(smin < s01 - 0.5*range) {
      smin = s01;
      below_min = true;
    }
    if(smax > s99 + 0.5*range) {
      smax = s99;
      above_max = true;
    }
    double denom = std::max(1.0, std::max(fabs(smax), fabs(smin)));
    if((smax-smin)/denom<1e-8) //min and max are almost the same
    smax = smin + 1e-8*denom;
    tmp.back() += 1e-8*denom;
    range = smax - smin;

    // plotting
    /*
    std::ofstream outfile;
    std::streambuf *buf;
    if(strcmp(iod_terminal.filename, "")) {
      if(first_time_open_file) {
        outfile.open(iod_terminal.filename, std::ofstream::out);
        first_time_open_file = false;
      } else
        outfile.open(iod_terminal.filename, std::ofstream::app);

      if(!outfile.is_open()) {
        fprintf(stdout,"\033[0;31m*** Error: Unable to open file %s for visualization.\n", iod_terminal.filename);
        exit(-1);
      }
      buf = outfile.rdbuf();
    } else
      buf = std::cout.rdbuf();

    std::ostream myout(buf);
*/
    *myout << "\n";
    *myout << "\033[0;32m==========================================\n";
    *myout << " " << solution_name << " at " << std::scientific << std::setprecision(4)
          << time << " (Step " << time_step << ")\n";
    *myout << "\033[0;32m==========================================\033[0m\n";

    if(iod_terminal.plane == TerminalVisualizationData::YZ)
      *myout << "The Y-Z plane. x: " << std::scientific << std::setprecision(2) << xyzmin[0] 
            << ", y (horiz.): [" << std::scientific << std::setprecision(2) << xyzmin[1] 
            << ", " << std::scientific << std::setprecision(2) << xyzmax[1]
            << "], z (vert.): [" << std::scientific << std::setprecision(2) << xyzmin[2] 
            << ", " << std::scientific << std::setprecision(2) << xyzmax[2] << "].\n";
    else if(iod_terminal.plane == TerminalVisualizationData::XZ)
      *myout << "The X-Z plane. y: " << std::scientific << std::setprecision(2) << xyzmin[1] 
            << ", x (horiz.): [" << std::scientific << std::setprecision(2) << xyzmin[0] 
            << ", " << std::scientific << std::setprecision(2) << xyzmax[0]
            << "], z (vert.): [" << std::scientific << std::setprecision(2) << xyzmin[2] 
            << ", " << std::scientific << std::setprecision(2) << xyzmax[2] << "].\n";
    else //XY
      *myout << "The X-Y plane. z: " << std::scientific << std::setprecision(2) << xyzmin[2] 
            << ", x (horiz.): [" << std::scientific << std::setprecision(2) << xyzmin[0] 
            << ", " << std::scientific << std::setprecision(2) << xyzmax[0]
            << "], y (vert.): [" << std::scientific << std::setprecision(2) << xyzmin[1] 
            << ", " << std::scientific << std::setprecision(2) << xyzmax[1] << "].\n";

    int color_code;
    double rel_value;
    std::vector<int>* ANSI_codes = (iod_terminal.colormap==TerminalVisualizationData::GRAYSCALE) ?
                                   &gray_ANSI : &turbo_ANSI;

    for(int i=0; i<ncols+1; i++) 
      *myout << "--";
    *myout << "\n";

    for(int j=0; j<nrows; j++) {
      *myout << "\u001b[0m|";
      for(int i=0; i<ncols; i++) {

        rel_value = (sol[j*ncols+i] - smin)/(smax - smin);
        if(rel_value<0)
          color_code = under_color_ANSI;
        else if(rel_value>1)
          color_code = over_color_ANSI;
        else 
          color_code = (*ANSI_codes)[round(rel_value*(ANSI_codes->size()-1))];
       
        *myout << "\u001b[48;5;" << color_code << "m  "; 
      }
      *myout << "\u001b[0m|\n";
    }
    for(int i=0; i<ncols+1; i++)
      *myout << "--";
    *myout << "\n";

    if(below_min)
      *myout << "\u001b[48;5;" << under_color_ANSI << "m " << std::scientific << std::setprecision(2)
            << tmp.front() << " ";
    int counter = 0;
    double step_size = (smax - smin)/(ANSI_codes->size()-1);
    double my_tick = smin;
    for(auto&& mycode : *ANSI_codes) {
      *myout << "\u001b[48;5;" << mycode << "m " << std::scientific << std::setprecision(2)
            << my_tick << " ";
      my_tick += step_size;
      if(++counter>8) {
        *myout << "\u001b[0m\n";
        counter = 0;
      }
    } 
    if(above_max)
      *myout << "\u001b[48;5;" << over_color_ANSI << "m " << std::scientific << std::setprecision(2)
            << tmp.back() << " ";

    *myout << "\u001b[0m\n" << std::endl;

    *myout << std::flush;

    if(outfile.is_open())
      outfile.flush();
  }

  MPI_Barrier(comm);

  if(strcmp(iod_terminal.filename, "") == 0) {
    usleep(1000*1000*iod_terminal.pause); //usleep counts in micro-seconds
  }

  iFrame++;
  last_snapshot_time = time;
  last_snapshot_clocktime = walltime();
  MPI_Allreduce(MPI_IN_PLACE, &last_snapshot_clocktime, 1, MPI_DOUBLE, MPI_MAX, comm); 

}

//------------------------------------------------------------------

