//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// Class Full_Morphology, 2D and 3D
//
//
// Applies Full Morphology Method based on the paper by Magnani et al, 2000 and
// several other references on Mathematical Morphology. Adapted from Zabot et al 2024.
//
//
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|

#ifndef FULL_MORPHOLOGY_HPP
#define FULL_MORPHOLOGY_HPP


#include <fstream>
#include <iomanip>
#include <stdint.h>
using namespace std;

#include "edt.hpp"
#include "geometry.hpp"
#include "component_labeling.hpp"

#include "file_uti.hpp"
#include "true_time.hpp"
#include "meusTipos.hpp"
#include "dataFileReader.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../common/Array.h"
#include "../common/Domain.h"
#include "../analysis/distance.h"

typedef Array<int> IntArray; //int32, LBPM type
// typedef Array<int16_t> Int16Array; 
typedef Array<bool> BoolArray;

typedef vector< vector<int> > MTvvi;






// -----------------------------------------------------------------------------
// Class
class Full_Morphology{
  public:
  
    Full_Morphology( int, char *[] );
    
    int Ndiam (void) const { return _d.size(); }
    int diameter( MTci & );
  
    void calc( MTci & );
    
    
  public:
  
    //int  _nthreads;
  
    int _ny, _nx, _nz;                                  // Work dimensions (with reservoirs, if any)
    int dimy, dimx, dimz;                               // Original dimensions
    IntArray _mmorig;                                   // Original image, used for reference
    IntArray _mm;                                       // Work image, modified during the simulation
    
    int _NP;                                            // Porous pixels
    
    MTvi _d;                                            // Diameters

    bool _compressible;                                 // Compressibility, set as true for MICP
    
    
    bool _iX, _iY, _iZ;                                 // Direction of invasion
    bool _iP, _iM;                                      // Sense of invasion
    int _chamber_i, _chamber_o;                         // Reservoir regions
    bool _all_faces;                                    // Surrounds the image for MICP   
    bool _saveImg;
    
    string _outImgRoot, _outImgDir; 
    ofstream _dat, _ftime;
    
    
    True_Time _ttime;                                   // Medição do tempo
    
  
    // Voxel labels

    int _I,                                             // I -> Inlet fluid
        _O,                                             // O -> Outlet fluid
        _S;                                             // S -> Solid

    
    BoolArray _trapped;                                 // Trapped regions
  
  
    bool has_out_inlet;
  
  
    IntArray _matrix1, _matrix2;
  
    IntArray _mx_edt1, _mx_edt2;

    IntArray _final_map;
   
    MTvi _next, _tail, _rtable;
  
    IntArray _edt;                                      // EDT matrix
  };
  
  
  
  
  
  
  
  
  
  
  //------------------------------------------------------------------------------
  // Description:
  //   Constructor
  //   Reads config and image file.
  //   Initializes class parameters.
  // INPUT:
  //   config => config file name
  //------------------------------------------------------------------------------
  Full_Morphology::Full_Morphology( int argc, char *argv[] ){


    char LocalRankFilename[40];
    string filename;
    double Rcrit_new;
    filename=argv[1];
    Rcrit_new=0.f; 
    NULL_USE( Rcrit_new );

    // Reads the input database 
    auto db = std::make_shared<Database>(filename);
    auto domain_db = db->getDatabase("Domain");
    auto fm_db = db->getDatabase("FM");
    
    auto size = domain_db->getVector<int>("N");
    _nx = size[0];
    _ny = size[1];
    _nz = size[2];
    
    auto VoxelLabels = domain_db->getVector<int>("VoxelLabels");
    _S = VoxelLabels[0];
    _O = VoxelLabels[1];
    _I = VoxelLabels[2];
    
    auto READFILE = domain_db->getScalar<std::string>("Filename");
    MTcs mmfile(READFILE);
    
    _outImgRoot = fm_db->getScalar<std::string>("ImageRoot");
    _outImgDir = fm_db->getScalar<std::string>("ImageDir");
    _saveImg = fm_db->getScalar<bool>("SaveImage");
    
    _all_faces = true;        //both true for MICP
    _compressible = true;
    
    if (!_all_faces) {
        auto direction = fm_db->getVector<int>("Direction");
        _iX = false;
        _iY = false;
        _iZ = false;
        _iM = false;
        _iP = false;
        
        if (direction[0] != 0 && direction[1] == 0 && direction[2] == 0) {
            _iX = true;
            if (direction[0] > 0) {
                _iP = true;
            } else {
                _iM = true;
            }
        } else if (direction[1] != 0 && direction[0] == 0 && direction[2] == 0) {
            _iY = true;
            if (direction[1] > 0) {
                _iP = true;
            } else {
                _iM = true;
            }
        } else if (direction[2] != 0 && direction[0] == 0 && direction[1] == 0) {
            _iZ = true;
            if (direction[2] > 0) {
                _iP = true;
            } else {
                _iM = true;
            }
        } else {
            aborta("Unkonwn direction.");
        }
    } else {
        cout << "Surround selected. Direction is disconsidered. " << endl;
        _iX = false;
        _iY = false;
        _iZ = false;
        _iM = false;
        _iP = true;

        if(_nz > 1)
            _iZ = true;
        else if(_ny > 1)
            _iY = true;
        else if(_nx > 1)
            _iX = true;

    }
    
    auto diameters = fm_db->getVector<int>("Diameters");
    for (int dd = diameters[0]; dd <= diameters[1]; dd += diameters[2]) {
        _d.push_back(dd);
    }
       sort( _d.begin(), _d.end() );
    
       MTvi::iterator it = unique( _d.begin(), _d.end() );
       _d.resize( distance( _d.begin(),it ) );
     
       if( _d.size()==0 )
       aborta("It was impossible to create diameters array.");
     
      reverse( _d.begin(), _d.end() );
    
  
    _outImgDir += "/";
    mymkdir( _outImgDir );
  
    
   
    double tt_read = _ttime();
  
  
    // Adds reservoirs where needed
    int x0=0, y0=0, z0=0;
    int xM=_nx, yM=_ny, zM=_nz;
    dimx = _nx; dimy = _ny; dimz = _nz; 
  
  

    if(!_all_faces){
    if     ( _iX ){ _nx += 2;   x0++;   xM = _nx-1; }
    else if( _iY ){ _ny += 2;   y0++;   yM = _ny-1; }
    else if( _iZ ){ _nz += 2;   z0++;   zM = _nz-1; }
    // ---------------------------------------------------------------------------
    }
    else {
      if( _nx>1 ){ _nx+=2;  x0++;  xM++; }
      if( _ny>1 ){ _ny+=2;  y0++;  yM++; }
      if( _nz>1 ){ _nz+=2;  z0++;  zM++; }
    }
  


    int iaux=-1;
    if     ( _iX ){ iaux=_nx-1; }
    else if( _iY ){ iaux=_ny-1; }
    else if( _iZ ){ iaux=_nz-1; }
  
    if( _iP ){
      _chamber_i =    0;
      _chamber_o = iaux;
    }else if( _iM ){
      _chamber_i = iaux;
      _chamber_o =    0;
    }
    // ---------------------------------------------------------------------------
   
  
  
  
    // ---------------------------------------------------------------------------
      _mmorig.resize(_nx, _ny, _nz);
      _mm.resize(_nx, _ny, _nz);
      _edt.resize(_nx, _ny, _nz);
      _matrix1.resize(_nx, _ny, _nz);
      _matrix2.resize(_nx, _ny, _nz);
      _mx_edt1.resize(_nx, _ny, _nz);
      _mx_edt2.resize(_nx, _ny, _nz);
      _trapped.resize(_nx, _ny, _nz);
      _final_map.resize(_nx, _ny, _nz);

    // MTci N3 = _nx*_ny*_nz;
    // MTci max_eq = static_cast<int> ( ceil( static_cast<double>(N3)/2.0 ) );       OLD

    size_t N3 = static_cast<size_t>(_nx) * static_cast<size_t>(_ny) * static_cast<size_t>(_nz);
    size_t max_eq = static_cast<size_t>(ceil(static_cast<double>(N3) / 2.0));

    _tail.resize(max_eq);
    _next.resize( max_eq );
    _rtable.resize( max_eq );
    // ---------------------------------------------------------------------------
   
  if(_all_faces) surround(_mm, _I);

    // ---------------------------------------------------------------------------
    // Inicializes matrices
    
    if(!_all_faces){
    // x invasion
    if( _iX ){
      for( int z=0; z<_nz; z++ ){
      for( int y=0; y<_ny; y++ ){
        _trapped(_chamber_i,y,z) = false;
        _trapped(_chamber_o,y,z) = false;
        
        _mm(_chamber_i,y,z) =_I;
        _mm(_chamber_o,y,z) =_O;
      }}
      
    // y invasion
    }else if( _iY ){
      for( int z=0; z<_nz; z++ ){
      for( int x=0; x<_nx; x++ ){
        _trapped(x,_chamber_i,z) = false;
        _trapped(x,_chamber_o,z) = false;
        
        _mm(x,_chamber_i,z) =_I;
        _mm(x,_chamber_o,z) =_O;
      }}
      
    // z invasion
    }else if( _iZ ){
      for( int y=0; y<_ny; y++ ){
      for( int x=0; x<_nx; x++ ){
        _trapped(x,y,_chamber_i) = false;
        _trapped(x,y,_chamber_o) = false;
        
        _mm(x,y,_chamber_i) = _I;
        _mm(x,y,_chamber_o) = _O;
      }}
      
    } }

    ifstream FRAW( mmfile.c_str() );
    abriu( FRAW, mmfile );
    unsigned char auxraw;
    int pc;

    has_out_inlet=false;
    for( int z=z0; z<zM; z++ ){
    for( int y=y0; y<yM; y++ ){
    for( int x=x0; x<xM; x++ ){
  
      _trapped(x,y,z) = false;
      FRAW >> auxraw;
      pc = static_cast<int>( auxraw );
      
      if( pc == _I )  has_out_inlet=true;
      if( pc != _S )  _NP++;
  
      if( pc!=_I  &&  pc!=_S  &&  pc!=_O ){  
        aborta("Unknown color in (" + ntos(x) + ", " +ntos(y) + ", " +ntos(z) + ")." );
      }
      
      _mm(x,y,z)     = pc;
    }
    }}
  
    // Copies _mm into _mmorig and initializes _final_map
    _NP=0;
    for( int z=0; z<_nz; z++ ){
    for( int y=0; y<_ny; y++ ){
    for( int x=0; x<_nx; x++ ){
      iaux = _mm(x,y,z);
      if( iaux!=_S ) _NP++;
      _mmorig(x,y,z) = iaux;

      if( iaux == _S ){
        _final_map(x,y,z) = 0;
      } else {
          _final_map(x,y,z) = -1;
        }

    }}}  
  
    int discount = _nx * _ny * _nz - dimx * dimy * dimz;

    // Discounts reservoir voxels

    if(!_all_faces){
    if     ( _iX ){ _NP -= 2*_ny*_nz;  }
    else if( _iY ){ _NP -= 2*_nx*_nz;  }
    else if( _iZ ){ _NP -= 2*_nx*_ny;  } }
    else{
      _NP = _NP - discount;
    }
  
  
  
    
    tt_read = _ttime() - tt_read;
    // ---------------------------------------------------------------------------
  
   
    
    // ---------------------------------------------------------------------------
    // Calculates _mmorig EDT
    double tt_edt = _ttime();
    euclidian_distance_transform( _S, _mmorig, _edt, _mx_edt1, _mx_edt2 );
    tt_edt = _ttime() - tt_edt;
    // ---------------------------------------------------------------------------
  
  
  
    // ---------------------------------------------------------------------------
    //Creates output file .dat
    MTvs vec(6), cmt(4);
    
    vec[0] = "Step";
    vec[1] = "Diameter (px)";
    vec[2] = "Number of pixels occupied by inlet fluid.";
    vec[3] = "Number of pixels occupied by inlet fluid / Number of porous pixels";
    vec[4] = "Number of pixels occupied by outlet fluid.";
    vec[5] = "Number of pixels occupied by outlet fluid  / Number of porous pixels";
    
  
    cmt[0] = "Image width  (x) (px)     : " + ntos( dimx );
    cmt[1] = "Image height (y) (px)     : " + ntos( dimy );
    cmt[2] = "Image planes (z) (px)     : " + ntos( dimz );
    cmt[3] = "Number of porous pixels   : " + ntos( _NP  ); 
  
    string saux = _outImgDir + "/" + _outImgRoot + ".dat";
    _dat.open( saux.c_str() );
    abriu( _dat, saux );
    outputFileHead( argc, argv, _dat, vec, cmt );  
    // ---------------------------------------------------------------------------
  
  
  
    // ---------------------------------------------------------------------------
    // Time file
    saux = _outImgDir + "/" + _outImgRoot + "_time.txt";
    _ftime.open( saux.c_str() );
    _ftime << scientific << setprecision(6);
    
    abriu( _ftime, saux );
    _ftime << "# File with execution times" << endl;
    _ftime << "# Times in seconds" << endl;
    _ftime << endl;
    
    _ftime << "# " << tt_read << " # Time to read and initialize data" << endl;
    _ftime << "# " << tt_edt  << " # Time to calculate EDT at initialization" << endl;
    
    
    _ftime << "\n\n"
           << "# Col  1: Step number\n"
           << "# Col  2: Step Diameter\n"
           << "# Col  3: tt_step\n"
           << "# Col  4: tt_loop1\n"
           << "# Col  5: tt_edt1\n"
           << "# Col  6: tt_loop2\n"
           << "# Col  7: tt_cl1\n"
           << "# Col  8: tt_loop3\n"
           << "# Col 9: tt_loop4\n"
           << "# Col 10: tt_loop5\n"
           << "# Col 11: tt_loop1c\n"
           << "# Col 12: tt_cl1c\n"
           << "# Col 13: tt_loop2c\n"
           << "# Col 14: tt_end\n"
           << endl;
    // ---------------------------------------------------------------------------  
  }
  
  
  
  
  
  
  
  
  
  
  
  //------------------------------------------------------------------------------
  // DESCRIPTION:
  //   Returns the diameter for a step
  // INPUT:
  //   step => Step number
  //------------------------------------------------------------------------------
  int Full_Morphology::diameter( MTci &step ){
    if( step<0 || step>=_d.size() )
      aborta( "Invalid step value." );
    return _d[step];
  }
  
  
  
  
  
  
  
  
  
  //------------------------------------------------------------------------------
  // DESCRIPTION:
  //   Calculates invasion for a given diameter
  // INPUT:
  //   step => current step number
  //------------------------------------------------------------------------------
  void Full_Morphology::calc( MTci &step ){
    int iaux;
    
    
    MTci D = this->diameter(step);
    _ftime << setw(3) << step << setw(4) << D; 
    double tt_step = _ttime();  
    MTcd D24 = D*D/4.0; 
  
    // Background and Foreground values for binary images
    MTci B=0, F=1;
  
  
    //----------------------------------------------------------------------------
    // Performs morphological opening operation to obtain H region from Magnani (eq. 5)
    //
    // The algorithm is based on the following Euclidean Distance Transform (EDT):
    // It is then straightforward to perform an erosion with a disc of radius
    // r simply by removing all pixels whose distance label is less than r. A
    // dilation is similarly performed by eroding the background.
    // described on the book IMAGE ANALYSIS FOR THE BIOLOGICAL SCIENCES, by
    // C A GLASBEY and G W HORGAN, chapter 5, page 10.
    //----------------------------------------------------------------------------
  
  
    //----------------------------------------------------------------------------
    // LOOP 1:
    // * Determines H region by opening:
    // - 1st erosion -> matrix1
    // * equals _mm to _mmorig
    //----------------------------------------------------------------------------
    double tt_loop1 = _ttime();
    // #pragma omp parallel for num_threads (_nthreads)
    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){

      _matrix2(x,y,z) = F;
      if( _edt(x,y,z)>=D24 ){
        _matrix1(x,y,z) = F;
        
      
      }else{
         _matrix1(x,y,z) = B;
        if(_mmorig(x,y,z)==_I ){
           _matrix1(x,y,z) = F;
           _matrix2(x,y,z) = B;  
        }
      }      
    }}}

    tt_loop1 = _ttime() - tt_loop1;
    //----------------------------------------------------------------------------
  
  
  
  
  
    //----------------------------------------------------------------------------
    // Find the connected regions
    int xaux=_nx/2, yaux=_ny/2, zaux=_nz/2;
    if     ( _iX ){ xaux=_chamber_i; }
    else if( _iY ){ yaux=_chamber_i; }
    else if( _iZ ){ zaux=_chamber_i; }
    int chamber_label=0;
    component_labeling( _matrix1, F, B, _next, _tail, _rtable ); 
    chamber_label =  _matrix1(xaux,yaux,zaux);
    // #pragma omp parallel for num_threads (_nthreads)
    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){
      if( _matrix1(x,y,z) == chamber_label ){
  
        _matrix1(x,y,z) = _matrix2(x,y,z);
        
      }else{
        _matrix1(x,y,z) = B;
      }
    }}}  
  
    //----------------------------------------------------------------------------
    // Now, compute the EDT of the negative, treating the Foreground pixels as
    // Background. Compute the EDT of the image _matrix1 and store the result in _matrix2.
    double tt_edt1 = _ttime();
    euclidian_distance_transform( F, _matrix1, _matrix2, _mx_edt1, _mx_edt2  );
    tt_edt1 = _ttime() - tt_edt1;
    //----------------------------------------------------------------------------
  

    //----------------------------------------------------------------------------
    // LOOP 2a:
    // * Determine H region by opening:
    // - 2nd erosion
    // * Highlights regions occupied by fluids in _mm
    //----------------------------------------------------------------------------
    double tt_loop2 = _ttime();
    // #pragma omp parallel for num_threads (_nthreads)
    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){
      if( _matrix2(x,y,z) < D24 ){
        _mm(x,y,z) = _I;

        if( _final_map(x,y,z) == -1 ){
          _final_map(x,y,z) = D;
        }

      }

      _matrix1(x,y,z) = (_mm(x,y,z)==_I)? F:B;
    }}}
    tt_loop2 = _ttime() - tt_loop2;
    
    
    
    //----------------------------------------------------------------------------
    // LOOP 2b:
    // * Disconnected regions of invaded fluid on matrix1
    //----------------------------------------------------------------------------
    double tt_cl1 = _ttime();
    component_labeling( _matrix1, F, B, _next, _tail, _rtable );
    tt_cl1 = _ttime() - tt_cl1;
    chamber_label =  _matrix1(xaux,yaux,zaux);  


    //----------------------------------------------------------------------------
    // LOOP 3:
    // * Replaces equivalent indices
    // * Determines G region
    // * Determines Omega region
    // * Determines what kind of fluid will be at each voxel according to Omega
    // * Set final image colors
    //----------------------------------------------------------------------------
      double tt_loop3 = _ttime();
    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){
  
      // G region: L U H   eq.6
      bool rG=false;
      if     ( _iX ){ rG  = (_mm(x,y,z)==_I || x==_chamber_i); }
      else if( _iY ){ rG  = (_mm(x,y,z)==_I || y==_chamber_i); }
      else if( _iZ ){ rG  = (_mm(x,y,z)==_I || z==_chamber_i); }
  
  
      // Operador K, generates Omega region   eq.11
      bool rO=false;
      if( rG && _matrix1(x,y,z) == chamber_label )
        rO = true;

      if( rO ){ 
        _mm(x,y,z) = _I;
      }else{
        if( _mm(x,y,z) != _S )
          _mm(x,y,z) = _O;
      }
      
    }}}
    tt_loop3 = _ttime() - tt_loop3;
    
  
  
    double tt_loop4 = _ttime();
    
    int x0=0  , y0=0  , z0=0;
    int xM=_nx, yM=_ny, zM=_nz;


    if(!_all_faces){
    if( _iX ){
      // #pragma omp parallel for num_threads (_nthreads)
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        _mm(_chamber_i,y,z) = _I;
        _mm(_chamber_o,y,z) = _O;
      }}
      
      x0=1;
      xM--;
         
    }else if( _iY ){
      // #pragma omp parallel for num_threads (_nthreads)
      for( int x=0; x<_nx; x++ ){
      for( int z=0; z<_nz; z++ ){
        _mm(x,_chamber_i,z) = _I;
        _mm(x,_chamber_o,z) = _O;
      }}
      
      y0=1;
      yM--;
      
    }else if( _iZ ){
      // #pragma omp parallel for num_threads (_nthreads)
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
        _mm(x,y,_chamber_i) = _I;
        _mm(x,y,_chamber_o)= _O;
      }}
      
      z0=1;
      zM--;    
    }
  }
    tt_loop4 = _ttime() - tt_loop4;      
     
    double tt_loop5 = _ttime();
    if( has_out_inlet ){
      // #pragma omp parallel for num_threads (_nthreads)
      for( int x=x0; x<xM; x++ ){
      for( int y=y0; y<yM; y++ ){
      for( int z=z0; z<zM; z++ ){
        if( _mmorig(x,y,z)==_I )
          _mm(x,y,z) = _I;
      }}}
    }
    tt_loop5 = _ttime() - tt_loop5;      
  
  
  
    //----------------------------------------------------------------------------
    // If the fluid is incompressible, determine the trapped regions:
    // 1 - Find expelled fluid regions disconnected to the expelled fluid reservoir
    // 2 - Set this and past trapped voxels
    // 3 - Save this pixels as expelled fluid
    //----------------------------------------------------------------------------
    double tt_loop1c=0;
    double tt_cl1c=0;
    double tt_loop2c=0;
    
    if( !_compressible ){
  
      //--------------------------------------------------------------------------
      // LOOP c1a:
        // * Determines if a pixel is trapped or not
      tt_loop1c = _ttime();
      // #pragma omp parallel for num_threads (_nthreads)
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){

        if( _trapped(x,y,z) )
          _mm(x,y,z) = _O;
  
        _matrix1(x,y,z) = (_mm(x,y,z)==_O)? F:B;
      }}}
      tt_loop1c = _ttime() - tt_loop1c;      
      
      tt_cl1c = _ttime();
      component_labeling( _matrix1, F, B, _next, _tail, _rtable );
      tt_cl1c = _ttime() - tt_cl1c;
      
    
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
  
      tt_loop2c = _ttime();
      // #pragma omp parallel for num_threads (_nthreads)
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){    
  
        if( _mm(x,y,z)==_O &&  _matrix1(x,y,z) != chamber_label )
          _trapped(x,y,z)=true;
      }}}
      tt_loop2c = _ttime() - tt_loop2c;
    }
  
        
    // ---------------------------------------------------------------------------
    // Creates .raw file
  
    bool createRAW=false;

    if( _saveImg==true && step==_d.size()-1 )
      createRAW=true;
    
    
    double tt_end = _ttime();

    uint8_t aux_raw;
    FILE *FRAW;
    if( createRAW ){
      
      xaux=_nx, yaux=_ny, zaux=_nz;
      if     ( _iX ){ xaux = _nx-2; }
      else if( _iY ){ yaux = _ny-2; }
      else if( _iZ ){ zaux = _nz-2; }
      
      int ndig = max(_d[0],_d[_d.size()-1]);
      string saux = ".invasion_diameters";
      MTcs fraw = _outImgDir + "/" + _outImgRoot + saux + ".raw";
  
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
      if     ( iaux==_I ) Ninlet++;
      else if( iaux==_O ) Noutlet++;
      
      if( createRAW ){
        int val = _final_map(x,y,z);
        int16_t val_save = static_cast<int16_t>(val);
        fwrite(&val_save, sizeof(int16_t), 1, FRAW);
      }
    }
    }
    }
    if( createRAW ) fclose(FRAW); 
    tt_end = _ttime() - tt_end;
    

    _dat << setprecision(6) << step << " " << D << " " << Ninlet << " " << Ninlet/(1.0*_NP) << " " << Noutlet << " " << Noutlet/(1.0*_NP) << endl;

    // Saves times
    tt_step = _ttime() - tt_step;  
    _ftime << setw(14) << tt_step
           << setw(14) << tt_loop1
           << setw(14) << tt_edt1
           << setw(14) << tt_loop2
           << setw(14) << tt_cl1
           << setw(14) << tt_loop3
           << setw(14) << tt_loop4
           << setw(14) << tt_loop5
           << setw(14) << tt_loop1c
           << setw(14) << tt_cl1c
           << setw(14) << tt_loop2c
           << setw(14) << tt_end
           << endl;
  }
  

  
#endif // FULL_MORPHOLOGY_HPP
