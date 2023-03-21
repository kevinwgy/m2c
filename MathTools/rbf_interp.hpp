#include <Eigen/Dense>

using namespace std;

namespace MathTools {

double vec_dot_product(int n, double a1[], double a2[]);

// basis functions
void phi1 ( int n, double r[], double r0, double v[] );
void phi2 ( int n, double r[], double r0, double v[] );
void phi3 ( int n, double r[], double r0, double v[] );
void phi4 ( int n, double r[], double r0, double v[] );

// interpolation functions
void rbf_interp ( int m, int nd, double xd[], double r0,
  void phi ( int n, double r[], double r0, double v[] ), double w[],
  int ni, double xi[], double fi[] );
void rbf_weight ( int m, int nd, double xd[], double r0,
  void phi ( int n, double r[], double r0, double v[] ),
  double fd[], double w[] );

} // end of namespace
