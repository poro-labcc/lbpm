//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// OBJETIVO:
//   Aplica o Young Laplace Method a uma geometria 2D ou 3D.
//________________________________________________________
// RECEBE:
//   nth => Número de threads para usar
//   cfg => Arquivo de configuração
//________________________________________________________
// MODO DE USAR:
//   ylm 4 config.txt
//________________________________________________________
//A.Z. - 12/13 => Criacao
//       02/14 => Adaptações para paralelização
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|


#include <iomanip>
#include <iostream>
using namespace std;

#include "../ylm/aborta.hpp"
#include "../ylm/numeros.hpp"
#include "../ylm/meusTipos.hpp"
#include "../ylm/ylm-comp.hpp"

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

  if( argc!=2 ) aborta("Wrong number of parameters.");
  Full_Morphology fm( argc, argv );
  
  MTci nsteps = fm.Ndiam();
  for( int step=0; step<nsteps; step++ ){
    
    MTci d = fm.diameter(step);
    cout << "Step " << step << ", D = " << d << " px." << endl;
    
    fm.calc(step);
    
  }

  return 0;
}











