//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
//
// Funções para ajudar na manipulação de arquivos
// Functions to help manipulate files
//________________________________________________________
//A.Z. - 03/05 => Creation
//       10/05 => Uses  aborta.h and meusTipos.h libraries
//       01/07 => Function outputFileHead
//       03/08 => Function outputFileHead also prints date
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef FILE_UTI_H
#define FILE_UTI_H

#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "aborta.hpp"
#include "meusTipos.hpp"
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
    aborta("It was impossible to open/create the file "+nome);
}




//------------------------------------------------------------------------------
// DESCRIPTION:
// Recieves a string that represents the name of a file and finds the path to the file.
//
// Atention: The string ist altered.
//
//
void substitui_caminho( string &arq ){

  int inicio  = arq.find("$");
  int fim     = arq.find("/");
  if(fim<0 && inicio>=0)
    fim = arq.size();
  int tamanho = fim - inicio - 1;

  if( inicio>=0 ){

    string variavel          = string( arq, inicio+1, tamanho );
    char *GETENV = getenv(variavel.c_str());

    if(GETENV==NULL){   
      string msg  = "Erro na funcao \"substitui_caminho(string &)\" da ";
      msg        += "biblioteca file_util.h\nNao foi possivel pegar a variavel";
      msg        += " de ambiente \"" + variavel + "\"";
      aborta(msg);
    }else              
      arq.replace(inicio,tamanho+1,string(GETENV));
  }
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
