//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//  Function that performs the Component Labeling form He, Chao & Suzuki, 2011.
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef COMPONENT_LABELING
#define COMPONENT_LABELING


#include <iostream>
using namespace std;


#include "fm_types.hpp"


#include "../common/Array.h"
typedef Array<int> IntArray;
typedef Array<bool> BoolArray;



void component_labeling( IntArray &IMG, MTci &F, MTci &B, MTvi &next,
                         MTvi &tail, MTvi &rtable );
void merge( MTci &, MTci &, MTvi &, MTvi &, MTvi & );
void resolve( MTci &, MTci &, MTvi &, MTvi &, MTvi & );






//------------------------------------------------------------------------------
// DESCRIPTION:
//   Identifies disconnected components of an image through He, Chao & Suzuki, 2011 algorithm.
//
// INPUTS:
//   img          => geometry to find disconnected regions
//   F, B         => Foreground and Background
//   next, tail, rtable => vector to work with equivalences information
//------------------------------------------------------------------------------
void component_labeling( IntArray &IMG, MTci &F, MTci &B, MTvi &next,
                         MTvi &tail, MTvi &rtable ){

  
  const int nx = IMG.size(0);
  const int ny = IMG.size(1);
  const int nz = IMG.size(2);
 
  int lx=0, nl=1;
  MTvi uniq_labels(3);
  int nuniq;
  
  
  for( int x=0; x<nx; x++ ){
  for( int y=0; y<ny; y++ ){
  for( int z=0; z<nz; z++ ){
    if( IMG(x,y,z)==F ){
      
      MTci lq = (x>0)?  IMG(x-1,y,z):B;
      MTci lp = (y>0)?  IMG(x,y-1,z):B;
      MTci lz = (z>0)?  IMG(x,y,z-1):B;

      nuniq=0;
      if( lp!=B ){
        uniq_labels[nuniq] = lp;
        nuniq++;
      }
      if( lq!=B && lq!=lp ){
        uniq_labels[nuniq] = lq;
        nuniq++;
      }
      if( lz!=B && lz!=lp && lz!=lq ){
        uniq_labels[nuniq] = lz;
        nuniq++;
      }


      switch( nuniq ){
    
        // ---------------------------------------------------------------------
        // Case 0: All neighbors are background
        case 0:
          nl++;
          lx = nl;
          
          // He, Chao, Suzuki, 2008, pág 752
          rtable[nl] = nl;
          next[nl]   = -1;
          tail[nl]   = nl;
          break;
      
      
        // ---------------------------------------------------------------------
        // Case 1: Only one label among the neighbors
        case 1:
          lx=uniq_labels[0];
          break;

      
        // ---------------------------------------------------------------------
        // Case 2: There are two different labels among the neighbors
        case 2:

          resolve( uniq_labels[0], uniq_labels[1], next, tail, rtable );
          
          lx = min( uniq_labels[0], uniq_labels[1] );
          break;

      
        // ---------------------------------------------------------------------
        // Case 3: There are three different labels among the neighbors
        case 3:

          resolve( uniq_labels[0], uniq_labels[1], next, tail, rtable );
          resolve( uniq_labels[0], uniq_labels[2], next, tail, rtable );
          resolve( uniq_labels[1], uniq_labels[2], next, tail, rtable );
          
          lx = min( uniq_labels[0], min( uniq_labels[1], uniq_labels[2] ) );
          break;
      }
      
      IMG(x,y,z) = lx;
      
            
    }
  }}}




  // Changes old indexes by equivalent ones
  for( int x=0; x<nx; x++ ){
  for( int y=0; y<ny; y++ ){
  for( int z=0; z<nz; z++ ){
    if( IMG(x,y,z)!=B )
      IMG(x,y,z) = rtable[ IMG(x,y,z) ];
  }}}


}



// He, Chao, Suzuki, 2008, page 752
void merge( MTci &u, MTci &v, MTvi &next, MTvi &tail, MTvi &rtable ){
  for( int i=v; i!=-1;  ){
    rtable[i] = u;
    i = next[i];
  }
  next[ tail[u] ] = v;
  tail[u] = tail[v];
}

// He, Chao, Suzuki, 2008, page 752
void resolve( MTci &x, MTci &y, MTvi &next, MTvi &tail, MTvi &rtable ){
  MTci u = rtable[x];
  MTci v = rtable[y];
  if     ( u<v ) merge( u, v, next, tail, rtable );
  else if( v<u ) merge( v, u, next, tail, rtable );
}


#endif // COMPONENT_LABELING