//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
//
//   Fornece abreviacoes uteis para os tipos mais usados
//   A regra de nomenclatura eh bem clara!
//
// PARA O EMACS:
//   Nao esqueca de alterar o arquivo emacs-lisp/font-lock.el para 
// incluir novos tipos. (Para achar o lugar, busque por MTci )
//   Para verificar se esta funcionando, basta fechar e reabrir o emacs com este
// arquivo. Se os novos tipos tiverem a mesma cor dos antigos, eh porque esta
// funcionando.
//________________________________________________________
//A.Z. - 07/05 => Criacao
//       12/05 => O comentario sobre o emacs.
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef MEUSTIPOS_H
#define MEUSTIPOS_H

#include <vector>
#include <string>
#include <iostream>
using namespace std;

typedef const int                   MTci;
typedef const bool                  MTcb;
typedef const long                  MTcl;
typedef const float                 MTcf;
typedef const double                MTcd;
typedef const string                MTcs;

typedef vector<int>                 MTvi;
typedef vector<bool>                MTvb;
typedef vector<long>                MTvl;
typedef vector<float>               MTvf;
typedef vector<double>              MTvd;
typedef vector<string>              MTvs;

typedef const vector<int>           MTcvi;
typedef const vector<bool>          MTcvb;
typedef const vector<long>          MTcvl;
typedef const vector<float>         MTcvf;
typedef const vector<double>        MTcvd;
typedef const vector<string>        MTcvs;

typedef vector<int>::iterator       MTvii;
typedef vector<bool>::iterator      MTvbi;
typedef vector<long>::iterator      MTvli;
typedef vector<float>::iterator     MTvfi;
typedef vector<double>::iterator    MTvdi;
typedef vector<string>::iterator    MTvsi;

#endif /* MEUSTIPOS_H */
