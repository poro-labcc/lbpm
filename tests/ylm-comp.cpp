//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//   Applies the Full Morphology Method to a 2D or 3D geometry.
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|


#include <iomanip>
#include <iostream>
using namespace std;

#include "../fm/fm_abort.hpp"
#include "../fm/fm_types.hpp"
#include "../fm/fm_method.hpp"

// #include "../fm/fm_numbers.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../common/Array.h"
#include "../common/Domain.h"
#include "../analysis/distance.h"
// #include "../analysis/morphology.h"









int main( int argc, char *argv[] ){

  if( argc!=2 ) abort_fm("Wrong number of parameters.");
  Full_Morphology fm( argc, argv );
  
  MTci nsteps = fm.Ndiam();
  for( int step=0; step<nsteps; step++ ){
    
    MTci d = fm.diameter(step);
    cout << "Step " << step << ", D = " << d << " px." << endl;
    
    fm.calc(step);
    
  }

  return 0;
}











