//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//  Funções que aplicam transformações na geometria das rochas.
//________________________________________________________
//A.Z. - 12/13 => Criação
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef FMM_GEOMETRY_HPP
#define FMM_GEOMETRY_HPP


// Padrão do C++
#include <iostream>
using namespace std;


// Minhas bibliotecas
#include "meusTipos.hpp"


// // Para matrizes 2D e 3D com Boost
// #include "boost/multi_array.hpp"
// typedef multi_array<int , 3> matrix;
// typedef multi_array<bool, 3> matrix_bool;

// // Desabilito as checagens do Boost para ganhar velocidade
// // Mas se houver algum problema, só vai dar pau no programa, nenhum aviso ;-)
// #define BOOST_DISABLE_ASSERTS

typedef Array<int> IntArray;
typedef Array<bool> BoolArray;



//------------------------------------------------------------------------------
// DESCRICAO:
//   Coloca paredes em uma geometria, conforme a direção de invasão
//
// RECEBE:
//   m => Matriz onde colocar as paredes
//   S => Número para colocar em um pixel sólido
//   a => Eixo
//
// ATENCAO:
//   m é alterada!
//
// RETORNA:
//   m => Matriz é alterada dentro da função, então retorna por referência
void walls( IntArray &m, MTci &S, MTcs &a ){
  const int nx = m.size(0);
  const int ny = m.size(1);
  const int nz = m.size(2);

  if( a=="x" ){ // Paredes em Y e Z
    if( ny > 1 ) for( int z=0; z<nz; z++ ) for( int x=0; x<nx; x++ ) { m(x,0,z) = S; m(x,ny-1,z) = S; }
    if( nz > 1 ) for( int y=0; y<ny; y++ ) for( int x=0; x<nx; x++ ) { m(x,y,0) = S; m(x,y,nz-1) = S; }
  } else if( a=="y" ){ // Paredes em X e Z
    if( nx > 1 ) for( int z=0; z<nz; z++ ) for( int y=0; y<ny; y++ ) { m(0,y,z) = S; m(nx-1,y,z) = S; }
    if( nz > 1 ) for( int y=0; y<ny; y++ ) for( int x=0; x<nx; x++ ) { m(x,y,0) = S; m(x,y,nz-1) = S; }
  } else if( a=="z" ){ // Paredes em X e Y
    if( nx > 1 ) for( int z=0; z<nz; z++ ) for( int y=0; y<ny; y++ ) { m(0,y,z) = S; m(nx-1,y,z) = S; }
    if( ny > 1 ) for( int z=0; z<nz; z++ ) for( int x=0; x<nx; x++ ) { m(x,0,z) = S; m(x,ny-1,z) = S; }
  }
}









//------------------------------------------------------------------------------
// DESCRICAO:
//   Coloca uma membrana semipermeável (1px) em uma geometria em um valor
//   específico de y, geralmente 0 ou ny-1.
//
// RECEBE:
//   m => Matriz onde colocar a membrana
//   S => Número para colocar em um pixel sólido
//   P => Número para colocar em um pixel poro
//   p => Posição onde colocar a membrana
//   a => Eixo
//
// ATENCAO:
//   m é alterada!
//
// RETORNA:
//   m => Matriz é alterada dentro da função, então retorna por referência
void membrane( IntArray &m, MTci &S, MTci &fill, MTci &p, MTcs &a ){
  const int nx = m.size(0);
  const int ny = m.size(1);
  const int nz = m.size(2);

  if ( a=="x" ){ // Membrana no plano X=p
    for( int z=0; z<nz; z++  ){
      for( int y=0; y<ny; y++ ) {
        if ((y+z)%2 == 0) m(p,y,z) = S;
        else m(p,y,z) = fill;
      }
    }
  } else if( a=="y" ){ // Membrana no plano Y=p
    for( int z=0; z<nz; z++  ){
      for( int x=0; x<nx; x++ ) {
        if ((x+z)%2 == 0) m(x,p,z) = S;
        else m(x,p,z) = fill;
      }
    }
  } else if( a=="z" ){ // Membrana no plano Z=p
    for( int y=0; y<ny; y++  ){
      for( int x=0; x<nx; x++ ) {
        if ((x+y)%2 == 0) m(x,y,p) = S;
        else m(x,y,p) = fill;
      }
    }
  }
}



void surround( IntArray &m, MTci &I){
  const int nx = m.size(0);
  const int ny = m.size(1);
  const int nz = m.size(2);

  // Preenche faces em X
  if (nx > 1) for(int z=0; z<nz; z++) for(int y=0; y<ny; y++) { m(0,y,z) = I; m(nx-1,y,z) = I; }
  // Preenche faces em Y
  if (ny > 1) for(int z=0; z<nz; z++) for(int x=0; x<nx; x++) { m(x,0,z) = I; m(x,ny-1,z) = I; }
  // Preenche faces em Z
  if (nz > 1) for(int y=0; y<ny; y++) for(int x=0; x<nx; x++) { m(x,y,0) = I; m(x,y,nz-1) = I; }
}



#endif // FMM_GEOMETRY_HPP





















