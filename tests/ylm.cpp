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

#include "../ylm/ylm.hpp"











int main( int argc, char *argv[] ){

  if( argc!=3 ) aborta("Número de parâmetros errado.");


  // Instancia e inicializa um objeto Full Morphology
  Young_Laplace_Method ylm( argc, argv );


  // Itera para cada diâmetro
  MTci nsteps = ylm.Ndiam();
  for( int step=0; step<nsteps; step++ ){
    
    MTci d = ylm.diameter(step);
    cout << "Passo " << step << ", D = " << d << " px." << endl;
    
    ylm.calc(step);
    
  }

  return 0;
}











