#include "Miniball_wrapper.h"
#include "Miniball.hpp"

using namespace foo;

void calc_miniball_wrapper(int dim, int n, double* points, double* center, double *radius){
  return foo::calc_miniball(dim, n, points, center, radius);
}
