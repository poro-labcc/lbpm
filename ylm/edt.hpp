//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//  Função que calcula a Transformada Euclidiana de Distâncias
//________________________________________________________
//A.Z. - 12/13 => Criação
//       02/14 => Parelização
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef EUCLIDIAN_DISTANCE_TRANSFORM
#define EUCLIDIAN_DISTANCE_TRANSFORM


// Padrão do C++
#include <iostream>
using namespace std;


// Minhas bibliotecas
#include "meusTipos.hpp"
//#include "true_time.hpp"


// Para matrizes 3D com Boost
#include "boost/multi_array.hpp"
using namespace boost;
typedef multi_array<int , 3> matrix;


// Para operações da Tranformada de Distância Euclidiana
#include "operators.hpp"










//------------------------------------------------------------------------------
// DESCRICAO:
//   Calcula Transformada de Distância Euclidiana 2D ou 3D
//   Aplico o algoritmo rápido (seção 3.4) de Saito & Toriwaki, 1994
//
//   http://www.sciencedirect.com/science/article/pii/0031320394901333
//
//   Baseado no código sedt.cc de David Coeurjolly
//   http://liris.cnrs.fr/~dcoeurjo
//
// ATENCAO:
//   O resultado é colocado numa matriz passada por referência por questão
//  de eficiência. Senão, eu retornaria, mas isso é lento.
//
//   Também por questão de eficiência, recebo duas matrizes auxiliares, ao
//  invés de criá-las a cada vez que a função é chamada.
//
// RECEBE:
//   background   => Cor do fundo
//   img          => Imagem para calcular a transformada
//   edt          => Matriz onde colocar o resultado da transformada
//   nthreads     => Número de threads para usar
//   maux, maux2  => Matrizes auxiliares já inicializadas com o tamanho da imagem
void euclidian_distance_transform( MTci &bg, const matrix &IMG, matrix &EDT,
                                   matrix &maux, matrix &maux2 ){
  
  //True_Time ttime;


  // Tamanho da imagem  
  MTci nx = IMG.shape()[0];
  MTci ny = IMG.shape()[1];
  MTci nz = IMG.shape()[2];


  
  // ---------------------------------------------------------------------------
  // First step of  the saito  algorithm
  //double tt_s1 = ttime();
  
  for( int y=0; y<ny; y++){
  for( int z=0; z<nz; z++){
    if( IMG[0][y][z] == bg ) maux[0][y][z] = 0;
    else                     maux[0][y][z] = INFTY;
  
    // Forward scan
    for( int x=1; x<nx; x++){
      if( IMG[x][y][z] == bg ) maux[x][y][z]= 0;
      else                     maux[x][y][z]= sum( 1, maux[x-1][y][z]);	  
    }
  
    //Backward scan
    for(int x=nx-2; x>=0; x--){    
      if( maux[x+1][y][z] < maux[x][y][z] ) 
        maux[x][y][z]=sum(1, maux[x+1][y][z]);
    }
  }}
  //tt_s1 = ttime() - tt_s1;





  // ---------------------------------------------------------------------------
  // Second step of the Saito algorithm using the
  //[Meijster/Roerdnik/Hesselink] optimization
  MTvi s(ny); //Center of the upper envelope parabolas
  MTvi t(ny); //Separating index between 2 upper envelope parabolas 
  int q, w;

  //double tt_s2 = ttime();

  for( int x=0; x<nx; x++){
  for( int z=0; z<nz; z++){
    q=0;
    s[0] = 0;
    t[0] = 0;

    //Forward Scan
    for( int u=1; u<ny; u++){
      
      while( (q >= 0) &&
       (F( t[q], s[q], prod(maux[x][s[q]][z],maux[x][s[q]][z]) ) > 
        F( t[q], u,    prod(maux[x][u][z],   maux[x][u][z])    ) )
           ){
        q--;
      }
    
      if(q<0){
        q=0;
        s[0]=u;
        
      }else{
        w = 1 + Sep( s[q],
                     u,
                     prod( maux[x][s[q]][z],maux[x][s[q]][z] ),
                     prod( maux[x][u][z]   ,maux[x][u][z]    ) );

        if( w < ny ){
          q++;
          s[q]=u;
          t[q]=w;
        }
      }
    }

    //Backward Scan
    for( int u=ny-1; u>=0; --u ){
      maux2[x][u][z] = F( u, s[q], prod( maux[x][s[q]][z], maux[x][s[q]][z]) );	   
      
      if( u==t[q] ) q--;
    }
  }}
  //tt_s2 = ttime() - tt_s2;


  // ---------------------------------------------------------------------------
  // Third step of the Saito algorithm using the
  //[Meijster/Roerdnik/Hesselink] optimization
  s.resize(nz); //Center of the upper envelope parabolas
  t.resize(nz); //Separating index between 2 upper envelope parabolas 

  //double tt_s3 = ttime();

  for( int x=0; x<nx; x++){
  for( int y=0; y<ny; y++){
    q=0;
    s[0] = 0;
    t[0] = 0;

    //Forward Scan
    for( int u=1; u<nz; u++){
      while( (q>=0) &&
             (F(t[q],s[q], maux2[x][y][s[q]]) > 
              F(t[q],u,maux2[x][y][u]))
           ){
             q--;
      }
      
      if(q<0){
        q=0;
        s[0]=u;
        
      }else{
        w = 1 + Sep( s[q],
                      u,
                      maux2[x][y][s[q]],
                      maux2[x][y][u] );
  
        if( w<nz ){
          q++;
          s[q]=u;
          t[q]=w;
        }
      }
    }
  
    //Backward Scan
    for( int u=nz-1; u>=0; --u){

      // Alterei o código aqui para ficar mais rápido
      // No original, havia a matriz sdt_xy
      // Mas não preciso dela, posso colocar direto no _edt        //sdt_xyz[x][y][u] = F( u, s[q], maux2[x][y][s[q]] );	      
      EDT[x][y][u] = F( u, s[q], maux2[x][y][s[q]] );	      

      if( u==t[q] ) 
        q--;
    }
  }}
  //tt_s3 = ttime() - tt_s3;


  //cout << "\n------------------------------------\n"
       //<< "Tempo 1: " << tt_s1 << "\n"
       //<< "Tempo 2: " << tt_s2 << "\n"
       //<< "Tempo 3: " << tt_s3 << "\n"
       //<< "Tempo T: " << tt_s1 + tt_s2 + tt_s3 << endl;
}







#endif // EUCLIDIAN_DISTANCE_TRANSFORM





















