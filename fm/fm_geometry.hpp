//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//  Functions that apply transformations to the geometries.
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|

#ifndef FMM_GEOMETRY_HPP
#define FMM_GEOMETRY_HPP



#include <iostream>
using namespace std;

#include "fm_types.hpp"


typedef Array<int> IntArray;   //LBPM TYPES
typedef Array<bool> BoolArray;



//------------------------------------------------------------------------------
// DESCRIPTION:
//   Sets walls in a geometry, according to invasion direction
//
// INPUTS:
//   m => matrix to be altered
//   S => voxel to be set at the wall
//   a => axis
//
//
// OUTPUT:
//   m => Matrix is altered inside the function, so it returns by reference
//------------------------------------------------------------------------------
void walls( IntArray &m, MTci &S, MTcs &a ){
  const int nx = m.size(0);
  const int ny = m.size(1);
  const int nz = m.size(2);

  if( a=="x" ){
    if( ny > 1 ) for( int z=0; z<nz; z++ ) for( int x=0; x<nx; x++ ) { m(x,0,z) = S; m(x,ny-1,z) = S; }
    if( nz > 1 ) for( int y=0; y<ny; y++ ) for( int x=0; x<nx; x++ ) { m(x,y,0) = S; m(x,y,nz-1) = S; }
  } else if( a=="y" ){
    if( nx > 1 ) for( int z=0; z<nz; z++ ) for( int y=0; y<ny; y++ ) { m(0,y,z) = S; m(nx-1,y,z) = S; }
    if( nz > 1 ) for( int y=0; y<ny; y++ ) for( int x=0; x<nx; x++ ) { m(x,y,0) = S; m(x,y,nz-1) = S; }
  } else if( a=="z" ){ 
    if( nx > 1 ) for( int z=0; z<nz; z++ ) for( int y=0; y<ny; y++ ) { m(0,y,z) = S; m(nx-1,y,z) = S; }
    if( ny > 1 ) for( int z=0; z<nz; z++ ) for( int x=0; x<nx; x++ ) { m(x,0,z) = S; m(x,ny-1,z) = S; }
  }
}









//------------------------------------------------------------------------------
// DESCRIPTION:
//   Sets a semipermeable membrane (1px) in a geometry and a specific coordinate
//   (usually 0 or _n -1)
//
// INPUTS:
//   m    => Matrix to be altered
//   S    => Solid voxel
//   fill => Fluid voxel 
//   p    => Membrane position
//   a    => Axis
//
// OUTPUT:
//   m => Matrix is altered inside the function, so it returns by reference
//------------------------------------------------------------------------------
void membrane( IntArray &m, MTci &S, MTci &fill, MTci &p, MTcs &a ){
  const int nx = m.size(0);
  const int ny = m.size(1);
  const int nz = m.size(2);

  if ( a=="x" ){
    for( int z=0; z<nz; z++  ){
      for( int y=0; y<ny; y++ ) {
        if ((y+z)%2 == 0) m(p,y,z) = S;
        else m(p,y,z) = fill;
      }
    }
  } else if( a=="y" ){
    for( int z=0; z<nz; z++  ){
      for( int x=0; x<nx; x++ ) {
        if ((x+z)%2 == 0) m(x,p,z) = S;
        else m(x,p,z) = fill;
      }
    }
  } else if( a=="z" ){
    for( int y=0; y<ny; y++  ){
      for( int x=0; x<nx; x++ ) {
        if ((x+y)%2 == 0) m(x,y,p) = S;
        else m(x,y,p) = fill;
      }
    }
  }
}


//------------------------------------------------------------------------------
// DESCRIPTION:
//   Adds a layer of voxels in all faces of the geometry (geometry must be surrounded)
//   by invaded fluid for MICP)
//
// INPUTS:
//   m    => Matrix to be altered
//   I    => Invaded fluid voxel
//
// OUTPUT:
//   m => Matrix is altered inside the function, so it returns by reference

void surround( IntArray &m, MTci &I){
  const int nx = m.size(0);
  const int ny = m.size(1);
  const int nz = m.size(2);

  if (nx > 1) for(int z=0; z<nz; z++) for(int y=0; y<ny; y++) { m(0,y,z) = I; m(nx-1,y,z) = I; }
  if (ny > 1) for(int z=0; z<nz; z++) for(int x=0; x<nx; x++) { m(x,0,z) = I; m(x,ny-1,z) = I; }
  if (nz > 1) for(int y=0; y<ny; y++) for(int x=0; x<nx; x++) { m(x,y,0) = I; m(x,y,nz-1) = I; }
}



#endif // FMM_GEOMETRY_HPP





















