//    Copright (C) 1999-2013, Bernd Gaertner
//    $Rev: 3581 $
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//    Contact:
//    --------
//    Bernd Gaertner
//    Institute of Theoretical Computer Science 
//    ETH Zuerich
//    CAB G31.1
//    CH-8092 Zuerich, Switzerland
//    http://www.inf.ethz.ch/personal/gaertner

// =========================================================================
// generates a set of random points, computes their smallest enclosing ball,
// and outputs the characteristics of this ball
// 
// usage: miniball_example [seed]
// =========================================================================

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "Miniball.hpp"



namespace foo {

  void calc_miniball(int dim, int n, double* points, double* center, double *radius){
  
    typedef double mytype;            // coordinate type  
  
    // for (int i=0;i<n*2;i++){
    //   printf("%f\n", p[i]);
    // }
  
    // generate random points and store them in a 2-d array
    // ----------------------------------------------------
    int itmp=0;
    mytype** ap = new mytype*[n];
    for (int i=0; i<n; ++i) {
      mytype* p = new mytype[dim];
      for (int j=0; j<dim; ++j) {
	p[j] = points[itmp]; //rand();

	//printf("%d,%f\n", itmp, points[itmp]);
	itmp++;      
      }
      ap[i]=p;
    }

    // define the types of iterators through the points and their coordinates
    // ----------------------------------------------------------------------
    typedef mytype* const* PointIterator; 
    typedef const mytype* CoordIterator;

    // create an instance of Miniball
    // ------------------------------
    typedef Miniball::
      Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > 
      MB;
    MB mb (dim, ap, ap+n);

    // output results
    // --------------
    // center
    const mytype* c = mb.center(); 
    for(int i=0; i<dim; ++i, ++c)
      center[i] = *c;

    // squared radius
    // std::cout << mb.squared_radius() <<  std::endl;
    *radius = std::sqrt(mb.squared_radius());
  
  }
}


  // int main (int argc, char* argv[])
  // {
  //   double p[] = {1.0, 0.0, 0.0, 1.0, -1.0, 0.0};
  //   double center[2];
  //   double radius;

  //   calc_miniball(2, 3, p, center, &radius);

  //   printf("%f,%f,%f\n", center[0], center[1], radius);
  // }

