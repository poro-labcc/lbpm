//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
//
// Funcao para abortar o programa quando ocorre um erro.
//________________________________________________________
//A.Z. - 03/05 => Criacao
//       11/05 => Usa #define ao inves de uma funcao.
//                Mostra o nome do arquivo e a linha onde o erro ocorreu
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef ABORTA_HPP
#define ABORTA_HPP


#include <iostream>
#include <cstdlib>
using namespace std;


#define aborta( erro ){                                         \
  cerr << "!!  Arquivo: " << __FILE__ << "\n"                   \
       << "!!  Linha  : " << __LINE__ << "\n"                   \
       << "!!  " << erro << "\n\a!!  ABORTANDO ..." << endl;    \
  exit(1);                                                      \
  }


#endif /* ABORTA_HPP */
