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

  // if( argc!=3 ) aborta("Número de parâmetros errado.");
  if( argc!=2 ) aborta("Wrong number of parameters.");

  // Instancia e inicializa um objeto Full Morphology
  Young_Laplace_Method ylm( argc, argv );

  // Initialize MPI
  // Utilities::startup( argc, argv );
  // Utilities::MPI comm( MPI_COMM_WORLD );
  //   int rank = comm.getRank();
  // {
  //   //.......................................................................
  //   // Reading the domain information file
  //   //.......................................................................
  //   char LocalRankFilename[40];
  //   string filename;
  //   double Rcrit_new;
  //   filename=argv[1];
  //   Rcrit_new=0.f; 
  //   NULL_USE( Rcrit_new );
  //   // read the input database 
	// 	auto db = std::make_shared<Database>( filename );
	// 	auto domain_db = db->getDatabase( "Domain" );

  //   auto size = domain_db->getVector<int>( "n" );
  //   ylm._nx = size[0];
  //   ylm._ny = size[1];
  //   ylm._nz = size[2];

  //   auto VoxelLabels = domain_db->getVector<int>( "VoxelLabels" );
  //   ylm._S = VoxelLabels[0];
  //   ylm._O = VoxelLabels[1];
  //   ylm._I = VoxelLabels[2];
  //   ylm._P = VoxelLabels[3];

  //   auto READFILE = domain_db->getScalar<std::string>( "Filename" );
  //   MTcs mmfile(READFILE);

  //   ylm._outImgRoot = domain_db->getScalar<std::string>( "ImageRoot" );
  //   ylm._outImgDir = domain_db->getScalar<std::string>( "ImageDir" );
  //   ylm._whichImg = domain_db->getScalar<std::string>( "WhichImage" );
  //   ylm._memb = domain_db->getScalar<int>( "Membrane" );
  //   ylm._wall = domain_db->getScalar<int>( "Walls" );
  //   ylm._wet = domain_db->getScalar<bool>( "Wetting" );
  //   ylm._compressible = domain_db->getScalar<bool>( "Compressible" );

  //   auto direction = domain_db->getScalar<char>( "Direction" );
  //   ylm._iX=false; ylm._iY=false; ylm._iZ=false;
  //   if(direction == 'x') {ylm._iX=true;}
  //   else if(direction == 'y') {ylm._iY=true;}
  //   else if(direction == 'z') {ylm._iZ=true;}
  //   else {aborta("Unkonwn direction.");}

  //   auto sense = domain_db->getScalar<char>( "Sense" );
  //   ylm._iM = false; ylm._iP = false;
  //   if(sense == '+'){ylm._iP = true;}
  //   else if(sense == '-'){ylm._iM = true;}
  //   else {aborta("Unknown sense.")}

  //   auto diameters = domain_db->getVector<int>( "Diameters" );
  //   for (int dd = diameters[0]; dd <= diameters[1]; dd+= diameters[2]){
  //     ylm._d.push_back(dd);
  //   }
  //      // Ordena
  //      sort( ylm._d.begin(), ylm._d.end() );
    
  //      // Tira repetidos
  //      MTvi::iterator it = unique( ylm._d.begin(), ylm._d.end() );
  //      ylm._d.resize( distance( ylm._d.begin(),it ) );
     
  //      // Aborta se não sobrou nada
  //      if( ylm._d.size()==0 )
  //      aborta("It was impossible to create diameters array.");
     
  //      // Se for não-molhante, reverte ordem do vetor
  //      if( !ylm._wet )  reverse( ylm._d.begin(), ylm._d.end() );



  // }

  // Itera para cada diâmetro
  MTci nsteps = ylm.Ndiam();
  for( int step=0; step<nsteps; step++ ){
    
    MTci d = ylm.diameter(step);
    cout << "Passo " << step << ", D = " << d << " px." << endl;
    
    ylm.calc(step);
    
  }

  return 0;
}











