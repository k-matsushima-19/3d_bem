#include "Miniball_wrapper.h"
#include <stdio.h>

int main(void){
  double p[] = {1.0, 0.0, 0.0, 1.0, -1.0, 0.0};
  double center[2];
  double radius;

  calc_miniball(2, 3, p, center, &radius);

  printf("%f,%f,%f\n", center[0], center[1], radius);
}
  
