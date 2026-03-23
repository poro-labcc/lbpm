//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//   Applies the Full Morphology Method to a 2D or 3D geometry.
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|


#include <iomanip>
#include <iostream>
using namespace std;

#include "../common/UtilityMacros.h"
#include "../fm/fm_method.hpp"

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

  if( argc!=2 ) ERROR("Wrong number of parameters.");
  Full_Morphology fm( argc, argv );
  
  const int nsteps = fm.Ndiam();
  for( int step=0; step<nsteps; step++ ){
    
    const int d = fm.diameter(step);
    cout << "Step " << step << ", D = " << d << " px." << endl;
    
    fm.calc(step);
    
  }

  return 0;
}










