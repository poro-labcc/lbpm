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
    int calc( const int & );
        
  public:

    std::vector<int> r_ini = {0,0,0};
    std::vector<int> r_end = {0,0,0};
    std::vector<int> iAxis = {false, false, false};     // If injection occour in the specified axis
    int _ny, _nx, _nz;                                  // Work dimensions (with reservoirs, if any)
    int dimy, dimx, dimz;                               // Original dimensions
    int _NP;                                            // Porous pixels
    int _chamber_i, _chamber_o;                         // Reservoir regions

    int rChamberI[3];
    int rChamberO[3];

    double resolution;
    
    vector<int> _diameter;                              // Diameter

    bool _compressible = false;                         // Compressibility, set as true for MICP
    bool _iP = false;                                   // Sense of invasion
    bool _all_faces;                                    // Surrounds the image for MICP   
    bool _saveImg;    
   
    UCharArray originalImage;                                // Original image, used for reference
    UCharArray workingImage;                                 // Work image, modified during the simulation   
    Array<int16_t> finalMap; 
    IntArray originalEDT;
    BoolArray trapped;                                               

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
    
    finalMap.resize( size[0],  size[1], size[2] );
    finalMap.fill(-1);

    auto ReadValues = domain_db->getVector<int>("ReadValues");
    auto WriteValues = domain_db->getVector<int>("WriteValues");

    resolution = domain_db->getScalar<double>("voxel_length");

    auto READFILE = domain_db->getScalar<std::string>("Filename");
    const string mmfile(READFILE);
    
    _saveImg = fm_db->getScalar<bool>("SaveImage");

    auto protocol = fm_db->getScalar<std::string>("protocol");
    
    checkOption( protocol, {"micp","drainage"} , "protocol" );
    if(protocol == "micp") _compressible=true; 
    _all_faces = _compressible;

    auto direction = fm_db->getScalar<std::string>("direction");
    checkOption( direction,  {"+x","-x","+y","-y","+z","-z","surround"} , "direction" );
    
    iAxis[0] = (direction[1] == 'x');
    iAxis[1] = (direction[1] == 'y');
    iAxis[2] = (direction[1] == 'z');

    _iP = (direction[0] == '+');

    if (direction == "surround") {
      _all_faces = true;
      _iP = true;
      if     (size[2] > 1) iAxis[2] = true;
      else if(size[1] > 1) iAxis[1] = true;
      else if(size[0] > 1) iAxis[0] = true; 
    }
    
    auto diameterRange = fm_db->getVector<int>("Diameters");
    int numberOfDiameters = (diameterRange[1] - diameterRange[0])/diameterRange[2] + 1;

    if(  numberOfDiameters <= 0 ) {
      ERROR("Error: It was impossible to create diameters array. ");  
    }

    _diameter.resize( numberOfDiameters  );

    for (int i = 0; i  < numberOfDiameters ; i++)  {
      _diameter[i] = diameterRange[1] - i * diameterRange[2];
    }

    r_end = size;
    dimx = _nx; dimy = _ny; dimz = _nz; 
    
    // Add aditional layer for input/output reservoirs    
    for  (int i = 0; i < 3; i++) {
         if ( ( (iAxis[i] ) && (!_all_faces) ) || ( (size[i]>1) && (_all_faces) ) ) {
            size[i] += 2;
            r_ini[i] = 1;
            r_end[i] = size[i] - 1;
         }          
    }

    _nx = size[0];
    _ny = size[1];
    _nz = size[2];
 
    originalImage.resize(_nx, _ny, _nz);
    workingImage.resize( originalImage.size() );
    originalEDT.resize( originalImage.size() );
    trapped.resize( originalImage.size() );

    workingImage.fill( INJECTED );  
    trapped.fill( false );
  
    // If injecting in a certain direction  then set the first or last faces as DISPLACED fluid
    for (int i = 0; i < 3; i++) {
      int rMin[3] = {0,0,0};
      int rMax[3] = {_nx,_ny,_nz};          

      if (iAxis[i] && (!_all_faces))
      {
          _chamber_i = _iP ? 0 :  size[i] -1 ;
          _chamber_o = _iP ? size[i] -1 : 0 ; 
          rMin[i] = _chamber_o;  
          rMax[i] = rMin[i] + 1;
          setRegion( workingImage     ,  DISPLACED  , rMin[0], rMax[0], rMin[1],rMax[1], rMin[2],rMax[2]);
      }
      
      rChamberI[i] = iAxis[i] ? _chamber_i : rMax[i] / 2;
      rChamberO[i] = iAxis[i] ? _chamber_o : rMax[i] / 2;
    }
    
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
    long expectedSize = (long) (r_end[2] - r_ini[2]) * (long) (r_end[1] - r_ini[1]) * (long) (r_end[0] - r_ini[0]);

    fseek(rawFile, 0, SEEK_END); 

    if ( ftell(rawFile) != expectedSize )
    {
        ERROR( "File '" +  mmfile + "' size is different from the expected (" + to_string(expectedSize) + " bytes)."  );
    }
    
    fseek(rawFile, 0, SEEK_BEGIN); // Move to beginning of the file to start reading
 
    unsigned char readValue;
    
    _NP= 0;
    for( int z= r_ini[2]; z< r_end[2]; z++ ) {
      for( int y = r_ini[1]; y< r_end[1]; y++ ) {
         for( int x= r_ini[0]; x< r_end[0]; x++ ) {

            fread(&readValue, sizeof( unsigned char), 1, rawFile);
            if (mapValue[readValue] == -1)             {
                ERROR( std::string("Not specified value in '" + filename + "' at (" +
                            to_string(x) + ", " + to_string(y) + ", " + to_string(z) + ")."  ) );
            }
        
            workingImage(x,y,z) = (unsigned char) mapValue[readValue];
            if( workingImage(x,y,z) != SOLID ) _NP++;
    }}}

    fclose(rawFile);
    originalImage = workingImage;
  
    // Calculates originalImage EDT
    edt_3d( SOLID, originalImage, originalEDT);

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
  
  int Full_Morphology::calc( const int &step ){
  
    const int D = _diameter[step];
    const double D24 = D*D/4.0;
       
    IntArray _matrix1 ( workingImage.size() );

    // Set as FOREGROUND the eroded POROUS space with diameter D united
    for( int z=0; z<_nz; z++ ) {
    for( int y=0; y<_ny; y++ ) {
    for( int x=0; x<_nx; x++ ) {

      if ( originalEDT(x,y,z) < D24 ) {              
        // Possible invaded region is selected as the one where the sphere fits and the one which it not located
        // at the DISPLACED fluid reservoir or SOLID 
        _matrix1(x,y,z) = ( originalImage(x,y,z) == INJECTED) ? FOREGROUND : BACKGROUND;
      }
      else _matrix1(x,y,z) = FOREGROUND;
    }}}
   
    // Label the eroded regions 
    component_labeling( _matrix1, FOREGROUND, BACKGROUND); 

    // Extract the label associated with the injection layer
    int chamber_label =  _matrix1( rChamberI[0],rChamberI[1],rChamberI[2]) ;

    for( int z=0; z<_nz; z++ ){
    for( int y=0; y<_ny; y++ ){
    for( int x=0; x<_nx; x++ ){
      if( _matrix1(x,y,z) == chamber_label ){
        _matrix1(x,y,z) = ( originalEDT(x,y,z) < D24 && originalImage(x,y,z) == INJECTED ) ? BACKGROUND : FOREGROUND;        
      }else{
        _matrix1(x,y,z) = BACKGROUND;
      }

    }}}  
  
    edt_3d( FOREGROUND, _matrix1, _matrix1);

    for( int z=0; z<_nz; z++ ){
    for( int y=0; y<_ny; y++ ){
    for( int x=0; x<_nx; x++ ){
    
      if( _matrix1(x,y,z) < D24 ){
        workingImage(x,y,z) = INJECTED;
      }
      _matrix1(x,y,z) = (workingImage(x,y,z)==INJECTED) ? FOREGROUND : BACKGROUND; 
    }}}
    
    component_labeling( _matrix1, FOREGROUND, BACKGROUND );
    chamber_label =  _matrix1( rChamberI[0], rChamberI[1], rChamberI[2]);  

    for( int z=0; z<_nz; z++ ){
    for( int y=0; y<_ny; y++ ){
    for( int x=0; x<_nx; x++ ){
        
        // G region: L U H   eq.6
      bool rG = (workingImage(x,y,z) == INJECTED);

      if     ( iAxis[0] ){ rG  = (rG || x == _chamber_i); }
      else if( iAxis[1] ){ rG  = (rG || y == _chamber_i); }
      else if( iAxis[2] ){ rG  = (rG || z == _chamber_i); }
    
      // Operador K, generates Omega region   eq.11
      bool rO = (rG && _matrix1(x,y,z) == chamber_label);
      
      if( rO ){ 
        workingImage(x,y,z) = INJECTED;
      }else{
        if( workingImage(x,y,z) != SOLID )
          workingImage(x,y,z) = DISPLACED;
      }
      
    }}}
  
    if(!_all_faces) {
      
      if ( iAxis[0] ){
          setRegion( workingImage     ,  INJECTED   , _chamber_i, _chamber_i+1, 0, _ny , 0, _nz );
          setRegion( workingImage     ,  DISPLACED  , _chamber_o, _chamber_o+1, 0, _ny , 0, _nz );      
      }
      else if( iAxis[1] ) {
          setRegion( workingImage     ,  INJECTED   , 0, _nx , _chamber_i, _chamber_i+1, 0, _nz );
          setRegion( workingImage     ,  DISPLACED  , 0, _nx , _chamber_o, _chamber_o+1, 0, _nz );           
      }
      else if( iAxis[2] ) {
          setRegion( workingImage     ,  INJECTED   , 0, _nx , 0, _nz , _chamber_i, _chamber_i+1);
          setRegion( workingImage     ,  DISPLACED  , 0, _nx , 0, _nz , _chamber_o, _chamber_o+1);                
      }
  }
         
    // If fluid is incompressible the disconect displaced regions must be keep
    // in the final state.
    if( !_compressible ){
  
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){

        if( trapped(x,y,z) )
          workingImage(x,y,z) = DISPLACED;
          _matrix1(x,y,z) = (workingImage(x,y,z)==DISPLACED)? FOREGROUND:BACKGROUND;

      }}}

      component_labeling( _matrix1, FOREGROUND, BACKGROUND );          
      chamber_label = _matrix1( rChamberO[0], rChamberO[1], rChamberO[2]);     

      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){    
        if( workingImage(x,y,z)==DISPLACED &&  _matrix1(x,y,z) != chamber_label ) trapped(x,y,z)=true;
      }}}
    
    }
  
    int injectedVolume = 0,  displacedVolume =0;
    for( int z= r_ini[2]; z< r_end[2]; z++ ) {
      for( int y = r_ini[1]; y< r_end[1]; y++ ) {
         for( int x= r_ini[0]; x< r_end[0]; x++ ) {

         if     ( workingImage(x,y,z) == INJECTED )  injectedVolume++;
         else if( workingImage(x,y,z) == DISPLACED ) displacedVolume++;

         int16_t* value = &finalMap(x - r_ini[0],y- r_ini[1],z- r_ini[2]);
         if( *value == -1 &&  workingImage(x,y,z) == INJECTED ) *value = (int16_t) D;
         if( *value == -1 && originalImage(x,y,z) == SOLID ) *value = 0;
    }

    }
    }

    if( _saveImg  && (D == _diameter.back() ) )
    {
      FILE *FRAW;  

      FRAW = fopen("invasion_diameters.raw", "wb");
      fwrite( finalMap.data() , sizeof(int16_t), finalMap.length(), FRAW);
      fclose(FRAW); 

      FILE *FMHD = fopen("invasion_diameters.mhd", "w");

      fprintf(FMHD, "ObjectType = Image\n");
      fprintf(FMHD, "NDims = 3\n");
      fprintf(FMHD, "DimSize = %d %d %d\n", dimx , dimy, dimz );
      fprintf(FMHD, "ElementType =  MET_SHORT\n");
      fprintf(FMHD, "ElementSpacing = %.1f %.1f %.1f\n", resolution, resolution, resolution);
      fprintf(FMHD, "ElementByteOrderMSB = False\n");
      fprintf(FMHD, "ElementDataFile = %s\n", "invasion_diameters.raw");
      fprintf(FMHD, "HeaderSize = 0\n");
      fclose(FMHD);
    }   

    FILE *log_file = fopen("injection_output.csv", "a");
    fprintf(log_file, "%d %d %f %d %f %d %f\n", step, D, D * resolution, injectedVolume, injectedVolume/(1.0*_NP),
          displacedVolume, displacedVolume/(1.0*_NP));
    fclose(log_file);

    return D;

  }
    
#endif // FULL_MORPHOLOGY_HPP
