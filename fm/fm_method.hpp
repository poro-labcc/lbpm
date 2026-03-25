//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// Class Full_Morphology, 2D and 3D
//
// Applies Full Morphology Method based on the paper by Magnani et al, 2000 and
// several other references on Mathematical Morphology. Adapted from Zabot et al 2024.
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|

#ifndef FULL_MORPHOLOGY_HPP
#define FULL_MORPHOLOGY_HPP

#define SOLID ((unsigned char) 0)
#define DISPLACED ((unsigned char) 1)
#define INJECTED ((unsigned char) 2)

#define BACKGROUND 0
#define FOREGROUND 1

#include <fstream>
#include <iomanip>
#include <stdint.h>
#include <string>
#include <bits/stdc++.h>

#include "fm_component_labeling.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../common/Array.h"
#include "../common/Domain.h"
#include "../analysis/distance.h"

using namespace std;

void checkOption( std::string a ,  std::vector<string> s, std::string keyName)
{
    std::string message = "Error: Invalid option '" + a + "' for " + keyName + ". Valid options are: " ;
    for (string b : s) {
            if (a == b) return;
            message += "'" + b + "', ";
    } 
    message.pop_back();
    ERROR( message );
}

template <typename TYPE>
void setRegion( Array<TYPE> &A, TYPE value, int x0, int x1, int y0, int y1, int z0, int z1)
{
    for (int z = z0; z < z1; z++) for (int y = y0; y < y1; y++) for (int x = x0; x < x1; x++) {
        A(x,y,z) = value;
    }   
}

static inline float intersection(int q, int vk, float fq, float fvk)
{
        float qq = (float)q  * (float)q;
        float vv = (float)vk * (float)vk;
        return ((fq + qq) - (fvk + vv)) / (2.0 * ((float)q - (float)vk));
}

void edt_1d(const int *f, int *g, int n)
{
    int *v = (int *) malloc((size_t)n * sizeof(int));
    float *z = (float *) malloc((size_t)(n + 1) * sizeof(float));

    int k = 0;
    v[0] = 0;
    z[0] = -INFINITY;
    z[1] =  INFINITY;

    for (int q = 1; q < n; q++) {
        float s;

        while (1) {
            int vk = v[k];
            s = intersection(q,vk,f[q],f[vk]);

            if (s <= z[k]) {
                k--;
                if (k < 0) {
                    k = 0;
                    break;
                }
            } else {
                break;
            }
        }

        if (k == 0) {
            int vk = v[k];
            s = intersection(q,vk,f[q],f[vk]);

            if (s <= z[k]) {
                v[0] = q;
                z[0] = -INFINITY;
                z[1] =  INFINITY;
                continue;
            }
        }

        k++;
        v[k] = q;
        z[k] = s;
        z[k + 1] = INFINITY;
    }

    k = 0;
    for (int x = 0; x < n; x++) {
        while (z[k + 1] < (float)x) {
            k++;
        }

        int vk = v[k];
        int dx = x - vk;
        g[x] = dx * dx + f[vk];
    }

    free(v);
    free(z);
}

void process_line(int* edt2, const size_t first, const size_t stride, int n)
{

            IntArray f(n);
            IntArray g(n);

            size_t pos = first;
            for (int i = 0; i < n; i++) {
                f(i) = edt2[pos];
                pos += stride;
            }

            edt_1d(f.data(), g.data(), n);

            pos = first;
            for (int i = 0; i < n; i++) {
                edt2[pos] = g(i);
                pos += stride;
            }
}

template <typename TYPE>
void edt_3d( unsigned char target, Array<TYPE> &image, IntArray &distance2)
{
    const int nx = image.size(0);
    const int ny = image.size(1);
    const int nz = image.size(2);

    size_t nvox =  image.length();
    int BIG = static_cast<float>(nx*nx + ny*ny + nz*nz) + 1;

    TYPE* img = image.data();
    int* edt2 = distance2.data();

    for (int i = 0; i < (int) nvox; i++) {
        edt2[i] = (img[i] == target) ? 0.0 : BIG;
    }

    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            process_line( edt2 , nx * (y + z * ny) , 1, nx );
        }
    }

    for (int z = 0; z < nz; z++) {
        for (int x = 0; x < nx; x++) {
            process_line( edt2 , z * nx * ny + x ,  nx, ny);
        }
    }

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            process_line( edt2 , y * nx + x , nx * ny, nz);
        }
    }
}

class Full_Morphology{
  public:
  
    Full_Morphology( int, char *[] );
    
    int Ndiam (void) const { return _d.size(); }
    int diameter( const int & );

    void calc( const int & );
        
  public:

    int _ny, _nx, _nz;                                  // Work dimensions (with reservoirs, if any)
    int dimy, dimx, dimz;                               // Original dimensions
    int _NP;                                            // Porous pixels
    int _chamber_i, _chamber_o;                         // Reservoir regions

    double resolution;
    
    vector<int> _d;                                            // Diameters

    bool _compressible = false;                                 // Compressibility, set as true for MICP
    bool _iX = false, _iY = false, _iZ = false;                                 // Direction of invasion
    bool _iP = false;                                   // Sense of invasion
    bool _all_faces;                                    // Surrounds the image for MICP   
    bool _saveImg;    
    bool has_out_inlet;

    string _outImgRoot;                         
    
    UCharArray _mmorig;                                   // Original image, used for reference
    UCharArray _mm;                                       // Work image, modified during the simulation
    IntArray _edt;
    BoolArray _trapped;                                               
  };

  Full_Morphology::Full_Morphology( int argc, char *argv[] ){

    string filename;
    
    filename=argv[1];  

    auto db = std::make_shared<Database>(filename);

    auto domain_db = db->getDatabase("Domain");
    auto fm_db = db->getDatabase("FM");
    
    auto size = domain_db->getVector<int>("N");
    _nx = size[0];
    _ny = size[1];
    _nz = size[2];
    
    auto ReadValues = domain_db->getVector<int>("ReadValues");
    auto WriteValues = domain_db->getVector<int>("WriteValues");

    resolution = domain_db->getScalar<double>("voxel_length");

    auto READFILE = domain_db->getScalar<std::string>("Filename");
    const string mmfile(READFILE);
    
    _outImgRoot = fm_db->getScalar<std::string>("ImageRoot");
    _saveImg = fm_db->getScalar<bool>("SaveImage");

    auto protocol = fm_db->getScalar<std::string>("protocol");
    
    checkOption( protocol, {"micp","drainage"} , "protocol" );
    if(protocol == "micp") _compressible=true; 
    _all_faces = _compressible;

    auto direction = fm_db->getScalar<std::string>("direction");
    checkOption( direction,  {"+x","-x","+y","-y","+z","-z","surround"} , "direction" );
    
    _iX = (direction[1] == 'x');
    _iY = (direction[1] == 'y');
    _iZ = (direction[1] == 'z');
    _iP = (direction[0] == '+');

    if (direction == "surround") {
      _all_faces = true;
      _iP = true;
      if(_nz > 1) _iZ = true;
      else if(_ny > 1) _iY = true;
      else if(_nx > 1) _iX = true; 
    }
    
    auto diameters = fm_db->getVector<int>("Diameters");

    int Ndiameters = (diameters[1] - diameters[0])/diameters[2] + 1;

    if(  Ndiameters <= 0 ) {
      ERROR("Error: It was impossible to create diameters array. ");  
    }

    _d.resize( Ndiameters );

    for (int i = 0; i  < Ndiameters; i++)  {
      _d[i] = diameters[1] - i * diameters[2];
    }

    std::vector<bool> iAxis = {_iX, _iY, _iZ };
    std::vector<int>  r_ini = {0,0,0};
    std::vector<int>  r_end = size;

    dimx = _nx; dimy = _ny; dimz = _nz; 
    
    // Add aditional layer for input/output reservoirs    
    for  (int i = 0; i < 3; i++) {
         if ( ( (iAxis[i] ) && (!_all_faces) ) || ( (size[i]>1) && (_all_faces) ) ) {
            size[i] += 2;
            r_ini[i] = 1;
            r_end[i] = size[i] - 1;
         }          
    }

    int x0 = r_ini[0] , y0 = r_ini[1]  , z0 =  r_ini[2];
    int xM = r_end[0],  yM = r_end[1]  , zM =  r_end[2];

    _nx = size[0];
    _ny = size[1];
    _nz = size[2];
 
    _mmorig.resize(_nx, _ny, _nz);
    _mm.resize(_nx, _ny, _nz);
    _mm.fill( INJECTED );

    _edt.resize(_nx, _ny, _nz);

    _trapped.resize(_nx, _ny, _nz);
    _trapped.fill( false );

  if (!_all_faces)  
  {     
    // If injecting in a certain direction  then set the first or last faces as DISPLACED fluid
    for (int i = 0; i < 3; i++)
    {
      int rMin[3] = {0,0,0};
      int rMax[3] = {_nx,_ny,_nz};
      if (iAxis[i])
      {
          _chamber_i = _iP ? 0 :  size[i] -1 ;
          _chamber_o = _iP ? size[i] -1 : 0 ; 
          rMin[i] = _chamber_o;  
          rMax[i] = rMin[i] + 1;
          setRegion( _mm     ,  DISPLACED  , rMin[0], rMax[0], rMin[1],rMax[1], rMin[2],rMax[2]);
      }
    }
   }

    has_out_inlet=false;
   
    int mapValue[255] = {-1};
    for(size_t idx = 0; idx < ReadValues.size(); idx++) {
      
          if ( (ReadValues[idx] < 0) || (ReadValues[idx] > 255) )
          {
            ERROR( "Only values between 0 - 255 can be used as labels in ReadValues");
            cout << ReadValues[idx] << endl;
          }
          if ( (WriteValues[idx] < 0) || (WriteValues[idx] > 2) )
          {
            ERROR( "Only values between 0 (SOLID), 1 and 2 (INJECT/DISPLACED FLUIDS) can be used as labels in WriteValues");
          }           
          mapValue[ ReadValues[idx] ] = (int) WriteValues[idx];           
    }

    FILE* rawFile = fopen( mmfile.c_str(), "r"); 
    if (rawFile == NULL) {
       ERROR("Error openning the file " + mmfile);
    }

    long SEEK_BEGIN = ftell(rawFile);
    long expectedSize = (long) (zM - z0) * (long) (yM - y0) * (long) (xM - x0);

    fseek(rawFile, 0, SEEK_END); 

    if ( ftell(rawFile) != expectedSize )
    {
        ERROR( "File '" +  mmfile + "' size is different from the expected (" + to_string(expectedSize) + " bytes)."  );
    }
    
    fseek(rawFile, 0, SEEK_BEGIN); // Move to beginning of the file to start reading
 
    unsigned char readValue;
    
    _NP= 0;
    for( int z=z0; z<zM; z++ ) {
      for( int y=y0; y<yM; y++ ) {
         for( int x=x0; x<xM; x++ ) {
            fread(&readValue, sizeof( unsigned char), 1, rawFile);
            if (mapValue[readValue] == -1)             {
                ERROR( std::string("Not specified value in '" + filename + "' at (" +
                            to_string(x) + ", " + to_string(y) + ", " + to_string(z) + ")."  ) );
            }
        
            _mm(x,y,z) = (unsigned char) mapValue[readValue];
            if( _mm(x,y,z) != SOLID ) _NP++;
    }}}

    fclose(rawFile);
    _mmorig = _mm;
  
    // Calculates _mmorig EDT
    edt_3d( SOLID, _mmorig, _edt);

    //Creates output file .csv and header
    bool WriteHeader = false;
    FILE *log_file = fopen("injection_output.csv", "r");
    if (log_file != NULL)
        fclose(log_file);
    else
        WriteHeader = true;

    if (WriteHeader) {
        log_file = fopen("injection_output.csv", "a+");
        fprintf(log_file, "step diameter_px diameter_um num_px_in frac_in num_px_out frac_out\n");
        fclose(log_file);
    }

  }
  
    //------------------------------------------------------------------------------
  // DESCRIPTION:
  //   Returns the diameter for a step
  // INPUT:
  //   step => Step number
  //------------------------------------------------------------------------------
  int Full_Morphology::diameter( const int &step ){
    if( step<0 || step>= (int) _d.size() )
      ERROR( "Invalid step value." );
    return _d[step];
  }
  
  //------------------------------------------------------------------------------
  // DESCRIPTION:
  //   Calculates invasion for a given diameter
  // INPUT:
  //   step => current step number
  //------------------------------------------------------------------------------
  void Full_Morphology::calc( const int &step ){
    int iaux;
  
    const int D = this->diameter(step);
    const double D24 = D*D/4.0;
       
    // Perform the opening of the pore space by a D diameter sphere, i. e.,
    // finds the region of porespace which can fit a D diameter sphere;
    IntArray _matrix1 ( _mm.size() );
    _matrix1.fill( FOREGROUND );

    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){      
      if ( _edt(x,y,z) < D24 ) {              
        // Possible invaded region is selected as the one where the sphere fits and the one which it not locate
        // at the DISPLACED fluid reservoir or SOLID 
        _matrix1(x,y,z) = ( _mmorig(x,y,z) != INJECTED) ? BACKGROUND : FOREGROUND;
      }
    }}}
   
    // Find the connected regions of the INJECTED FLUID
    component_labeling( _matrix1, FOREGROUND, BACKGROUND); 

    int xaux=_nx/2, yaux=_ny/2, zaux=_nz/2;
    if     ( _iX ){ xaux=_chamber_i; }
    else if( _iY ){ yaux=_chamber_i; }
    else if( _iZ ){ zaux=_chamber_i; }
    int chamber_label=0;

    chamber_label =  _matrix1(xaux,yaux,zaux);
    // #pragma omp parallel for num_threads (_nthreads)
    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){
      if( _matrix1(x,y,z) == chamber_label ){
        _matrix1(x,y,z) = ( _edt(x,y,z) < D24 && _mmorig(x,y,z) == INJECTED ) ? BACKGROUND : FOREGROUND;        
      }else{
        _matrix1(x,y,z) = BACKGROUND;
      }
    }}}  
  
    edt_3d( FOREGROUND, _matrix1, _matrix1);
    //----------------------------------------------------------------------------
  

    //----------------------------------------------------------------------------
    // LOOP 2a:
    // * Determine H region by opening:
    // - 2nd erosion
    // * Highlights regions occupied by fluids in _mm
    //----------------------------------------------------------------------------
    // #pragma omp parallel for num_threads (_nthreads)
    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){
      if( _matrix1(x,y,z) < D24 ){
        _mm(x,y,z) = INJECTED;
      }

      _matrix1(x,y,z) = (_mm(x,y,z)==INJECTED)? FOREGROUND : BACKGROUND; 
    }}}
    
    
    
    //----------------------------------------------------------------------------
    // LOOP 2b:
    // * Disconnected regions of invaded fluid on matrix1
    //----------------------------------------------------------------------------
    component_labeling( _matrix1, FOREGROUND, BACKGROUND );
    chamber_label =  _matrix1(xaux,yaux,zaux);  


    //----------------------------------------------------------------------------
    // LOOP 3:
    // * Replaces equivalent indices
    // * Determines G region
    // * Determines Omega region
    // * Determines what kind of fluid will be at each voxel according to Omega
    // * Set final image colors
    //----------------------------------------------------------------------------
    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){
  
      // G region: L U H   eq.6
      bool rG=false;
      if     ( _iX ){ rG  = (_mm(x,y,z)==INJECTED || x==_chamber_i); }
      else if( _iY ){ rG  = (_mm(x,y,z)==INJECTED || y==_chamber_i); }
      else if( _iZ ){ rG  = (_mm(x,y,z)==INJECTED || z==_chamber_i); }
  
  
      // Operador K, generates Omega region   eq.11
      bool rO=false;
      if( rG && _matrix1(x,y,z) == chamber_label )
        rO = true;

      if( rO ){ 
        _mm(x,y,z) = INJECTED;
      }else{
        if( _mm(x,y,z) != SOLID )
          _mm(x,y,z) = DISPLACED;
      }
      
    }}}
    
    int x0=0  , y0=0  , z0=0;
    int xM=_nx, yM=_ny, zM=_nz;


    if(!_all_faces){
    if( _iX ){
      // #pragma omp parallel for num_threads (_nthreads)
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        _mm(_chamber_i,y,z) = INJECTED;
        _mm(_chamber_o,y,z) = DISPLACED;
      }}
      
      x0=1;
      xM--;
         
    }else if( _iY ){
      // #pragma omp parallel for num_threads (_nthreads)
      for( int x=0; x<_nx; x++ ){
      for( int z=0; z<_nz; z++ ){
        _mm(x,_chamber_i,z) = INJECTED;
        _mm(x,_chamber_o,z) = DISPLACED;
      }}
      
      y0=1;
      yM--;
      
    }else if( _iZ ){
      // #pragma omp parallel for num_threads (_nthreads)
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
        _mm(x,y,_chamber_i) = INJECTED;
        _mm(x,y,_chamber_o)= DISPLACED;
      }}
      
      z0=1;
      zM--;    
    }
  }
    if( has_out_inlet ){
      // #pragma omp parallel for num_threads (_nthreads)
      for( int x=x0; x<xM; x++ ){
      for( int y=y0; y<yM; y++ ){
      for( int z=z0; z<zM; z++ ){
        if( _mmorig(x,y,z)==INJECTED )
          _mm(x,y,z) = INJECTED;
      }}}
    }
  
  
  
    //----------------------------------------------------------------------------
    // If the fluid is incompressible, determine the trapped regions:
    // 1 - Find expelled fluid regions disconnected to the expelled fluid reservoir
    // 2 - Set this and past trapped voxels
    // 3 - Save this pixels as expelled fluid
    //----------------------------------------------------------------------------

    if( !_compressible ){
  
      //--------------------------------------------------------------------------
      // LOOP c1a:
        // * Determines if a pixel is trapped or not
      // #pragma omp parallel for num_threads (_nthreads)
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){

        if( _trapped(x,y,z) )
          _mm(x,y,z) = DISPLACED;
  
        _matrix1(x,y,z) = (_mm(x,y,z)==DISPLACED)? FOREGROUND:BACKGROUND;
      }}}
      component_labeling( _matrix1, FOREGROUND, BACKGROUND );
      
    
      xaux=_nx/2, yaux=_ny/2, zaux=_nz/2;
      if     ( _iX ){ xaux=_chamber_o; }
      else if( _iY ){ yaux=_chamber_o; }
      else if( _iZ ){ zaux=_chamber_o; }
      chamber_label = _matrix1(xaux,yaux,zaux);
      
      

      //--------------------------------------------------------------------------
      // LOOP c2:
      // * Replaces equivalent indices
      // * All expelled fluid regions disconnected to the outlet are set as trapped.
      //--------------------------------------------------------------------------
  
      // #pragma omp parallel for num_threads (_nthreads)
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){    
  
        if( _mm(x,y,z)==DISPLACED &&  _matrix1(x,y,z) != chamber_label )
          _trapped(x,y,z)=true;
      }}}
    }
  
        
    // ---------------------------------------------------------------------------
    // Creates .raw file
  
    bool createRAW=false;

    if( _saveImg==true && step== (int) _d.size()-1 )
      createRAW=true;
    
    

    FILE *FRAW;
    if( createRAW ){
      
      xaux=_nx, yaux=_ny, zaux=_nz;
      if     ( _iX ){ xaux = _nx-2; }
      else if( _iY ){ yaux = _ny-2; }
      else if( _iZ ){ zaux = _nz-2; }
      
      string saux = "invasion_diameters";
      const string fraw = saux + ".raw";
  
      FRAW = fopen64(fraw.c_str(), "wb");
    }
        
        
    // ---------------------------------------------------------------------------
    // Last loop to save colors and count voxels occupied by each fluid
    x0=0  , y0=0  , z0=0;
    xM=_nx, yM=_ny, zM=_nz;

    if(!_all_faces){
    if     ( _iX ){ x0=1; xM = _nx-1; }
    else if( _iY ){ y0=1; yM = _ny-1; }
    else if( _iZ ){ z0=1; zM = _nz-1; } }
    else{
      if (dimx > 1){x0=1; xM = _nx-1;}
      if (dimy > 1){y0=1; yM = _ny-1;}
      if (dimz > 1){z0=1; zM = _nz-1;}
      }

    int Ninlet=0, Noutlet=0;
    for( int z=z0; z<zM; z++ ){
    for( int y=y0; y<yM; y++ ){
    for( int x=x0; x<xM; x++ ){    
      
      iaux = _mm(x,y,z);

      // if( _final_map(x,y,z) == -1 && _mm(x,y,z) == INJECTED ){
      //   _final_map(x,y,z) = D;
      // }

      if     ( iaux==INJECTED ) Ninlet++;
      else if( iaux==DISPLACED ) Noutlet++;
      
      // if( createRAW ){
      //   int val = _final_map(x,y,z);
      //   int16_t val_save = static_cast<int16_t>(val);
      //   fwrite(&val_save, sizeof(int16_t), 1, FRAW);
      // }
    }
    }
    }
    if( createRAW ) fclose(FRAW); 
    

//---------------------------------------------------------------
//fills in .csv
  FILE *log_file = fopen("injection_output.csv", "a");
  fprintf(log_file, "%d %d %f %d %f %d %f\n", step, D, D * resolution, Ninlet, Ninlet/(1.0*_NP),
         Noutlet, Noutlet/(1.0*_NP));
  fclose(log_file);


  }
  
  
#endif // FULL_MORPHOLOGY_HPP

