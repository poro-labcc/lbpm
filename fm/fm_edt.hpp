//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//  Function that calculates the EDT - Euclidian Distance Transform
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef EUCLIDIAN_DISTANCE_TRANSFORM
#define EUCLIDIAN_DISTANCE_TRANSFORM


#include <iostream>
using namespace std;


#include "fm_types.hpp"


#include "../common/Array.h"
typedef Array<int> IntArray;
typedef Array<bool> BoolArray;

#include "fm_operators.hpp"










//------------------------------------------------------------------------------
// DESCRIPTION:
//   Calculates 2D or 3D EDT
//   Applies algorithm (section 3.4) form Saito & Toriwaki, 1994
//
//   http://www.sciencedirect.com/science/article/pii/0031320394901333
//
//   Based on sedt.cc from David Coeurjolly
//   http://liris.cnrs.fr/~dcoeurjo
//
//
// INPUTS:
//   background   => background color
//   img          => image to calculate EDT
//   edt          => matrix to save EDT results
//   maux, maux2  => auxiliary matrices
void euclidian_distance_transform( MTci &bg, const IntArray &IMG, IntArray &EDT,
                                   IntArray &maux, IntArray &maux2 ){

    const int nx = IMG.size(0);
    const int ny = IMG.size(1);
    const int nz = IMG.size(2);
  
    // ---------------------------------------------------------------------------
    // STEP 1: X axis
    for( int y=0; y<ny; y++){
    for( int z=0; z<nz; z++){
        if( IMG(0,y,z) == bg ) maux(0,y,z) = 0;
        else                     maux(0,y,z) = INFTY;
      
        // forward scan
        for( int x=1; x<nx; x++){
            if( IMG(x,y,z) == bg ) maux(x,y,z) = 0;
            else                     maux(x,y,z) = sum( 1, maux(x-1,y,z));   
        }
      
        // backward scan
        for(int x=nx-2; x>=0; x--){    
            if( maux(x+1,y,z) < maux(x,y,z) ) 
                maux(x,y,z) = sum(1, maux(x+1,y,z));
        }
    }}

    // ---------------------------------------------------------------------------
    // STEP 2: Y axis
    MTvi s(ny > nz ? ny : nz); 
    MTvi t(ny > nz ? ny : nz); 
    int q, w;

    for( int x=0; x<nx; x++){
    for( int z=0; z<nz; z++){
        q=0;
        s[0] = 0;
        t[0] = 0;
    
        // forward scan
        for( int u=1; u<ny; u++){ 
            while( (q >= 0) &&
             (F( t[q], s[q], prod(maux(x,s[q],z),maux(x,s[q],z)) ) > 
              F( t[q], u,    prod(maux(x,u,z),   maux(x,u,z))    ) )
                 ){
                q--;
            }
        
            if(q<0){
                q=0;
                s[0]=u;
            } else {
                w = 1 + Sep( s[q],
                             u,
                             prod( maux(x,s[q],z), maux(x,s[q],z) ),
                             prod( maux(x,u,z),    maux(x,u,z)    ) );
    
                if( w < ny ){
                    q++;
                    s[q]=u;
                    t[q]=w;
                }
            }
        }
    
        // backward scan
        for( int u=ny-1; u>=0; --u ){
            maux2(x,u,z) = F( u, s[q], prod( maux(x,s[q],z), maux(x,s[q],z)) );    
            if( u==t[q] ) q--;
        }
    }}

    // ---------------------------------------------------------------------------
    // STEP 3: Z axis
    for( int x=0; x<nx; x++){
    for( int y=0; y<ny; y++){
        q=0;
        s[0] = 0;
        t[0] = 0;
    
        // forward scan
        for( int u=1; u<nz; u++){
            while( (q>=0) &&
                 (F(t[q],s[q], maux2(x,y,s[q])) > 
                  F(t[q],u,   maux2(x,y,u)))
                 ){
                 q--;
            }
          
            if(q<0){
                q=0;
                s[0]=u;
            } else {
                w = 1 + Sep( s[q],
                              u,
                              maux2(x,y,s[q]),
                              maux2(x,y,u) );
      
                if( w<nz ){
                    q++;
                    s[q]=u;
                    t[q]=w;
                }
            }
        }
      
        // backward scan
        for( int u=nz-1; u>=0; --u){
            EDT(x,y,u) = F( u, s[q], maux2(x,y,s[q]) );       
            if( u==t[q] ) 
                q--;
        }
    }}
}







#endif // EUCLIDIAN_DISTANCE_TRANSFORM





















