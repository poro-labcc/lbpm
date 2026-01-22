//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
//
// Functions to help manipulate files
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|

// THIS FILE IS NO MORE NECESSARY.

#ifndef FILE_UTI_H
#define FILE_UTI_H

#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "fm_abort.hpp"
#include "fm_types.hpp"
using namespace std;


void mymkdir( MTcs & );
template < class T > inline void abriu( T &, MTcs & );
void substitui_caminho( string & );
void outputFileHead( MTci &, char **, ofstream &, MTcvs & );




//------------------------------------------------------------------------------
// DESCRIPTION:
//   Creates a folder
// INPUT:
//   folder => folder name
void mymkdir( MTcs &folder ){
  MTcs mkdir( "mkdir -p " + folder );
  int lixo = system( mkdir.c_str() );
  lixo=lixo;
}




//------------------------------------------------------------------------------
// DESCRIPTION:
// Verifies if a file was open correctly
template < class T >
inline void abriu( T &obj, MTcs &nome ){
  if( !obj.is_open() )
    abort_fm("Unable to open/create the file "+nome);
}

//------------------------------------------------------------------------------
// DESCRIPTION:
//   Generates header on the output file.
void outputFileHead( MTci &argc, char *argv[], ofstream &OUT, MTcvs &vec,
                     MTcvs &comment){
  char *dir = getenv("PWD");

  OUT<< "# ##################################################################\n"
     << "# File created with the command:  ";
  for( int i=0;i<argc;i++ ) OUT << argv[i] << " ";

  OUT<< "\n# On folder:  ";
  if( dir != NULL ) OUT << dir;
  else              OUT << "(It was impossible to identify the folder!)";

  OUT<< "\n# On:  " << __TIME__ << "  " << __DATE__  << endl;


  OUT<< "#\n";
  if( !comment.empty() ){
    for( unsigned int i=0;i<comment.size();i++ )
      OUT << "#   * " << comment[i] << "\n";
    OUT<< "#\n";
  }

  for( unsigned int i=0;i<vec.size();i++ )
    OUT << "# Column " << setw(2) << i+1 << ": " << vec[i] << "\n";
  OUT<< "# ##################################################################"
     << endl;
}



#endif /* FILE_UTI_H */
