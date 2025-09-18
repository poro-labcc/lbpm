//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//  Função que calcula a Transformada Euclidiana de Distâncias
//________________________________________________________
//A.Z. - 12/13 => Criação
//       02/14 => Parelização
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef EUCLIDIAN_DISTANCE_TRANSFORM
#define EUCLIDIAN_DISTANCE_TRANSFORM


// Padrão do C++
#include <iostream>
using namespace std;


// Minhas bibliotecas
#include "meusTipos.hpp"
//#include "true_time.hpp"


// // Para matrizes 3D com Boost
// #include "boost/multi_array.hpp"
// using namespace boost;
// typedef multi_array<int , 3> matrix;

#include "../common/Array.h"
typedef Array<int> IntArray;
typedef Array<bool> BoolArray;

// Para operações da Tranformada de Distância Euclidiana
#include "operators.hpp"










//------------------------------------------------------------------------------
// DESCRICAO:
//   Calcula Transformada de Distância Euclidiana 2D ou 3D
//   Aplico o algoritmo rápido (seção 3.4) de Saito & Toriwaki, 1994
//
//   http://www.sciencedirect.com/science/article/pii/0031320394901333
//
//   Baseado no código sedt.cc de David Coeurjolly
//   http://liris.cnrs.fr/~dcoeurjo
//
// ATENCAO:
//   O resultado é colocado numa matriz passada por referência por questão
//  de eficiência. Senão, eu retornaria, mas isso é lento.
//
//   Também por questão de eficiência, recebo duas matrizes auxiliares, ao
//  invés de criá-las a cada vez que a função é chamada.
//
// RECEBE:
//   background   => Cor do fundo
//   img          => Imagem para calcular a transformada
//   edt          => Matriz onde colocar o resultado da transformada
//   nthreads     => Número de threads para usar
//   maux, maux2  => Matrizes auxiliares já inicializadas com o tamanho da imagem
// Esta versão foi traduzida para usar IntArray com acesso (x, y, z)
void euclidian_distance_transform( MTci &bg, const IntArray &IMG, IntArray &EDT,
                                   IntArray &maux, IntArray &maux2 ){
  
    // Lê as dimensões na ordem correta (x, y, z)
    const int nx = IMG.size(0);
    const int ny = IMG.size(1);
    const int nz = IMG.size(2);
  
    // ---------------------------------------------------------------------------
    // PASSO 1: Varredura ao longo do eixo X
    for( int y=0; y<ny; y++){
    for( int z=0; z<nz; z++){
        if( IMG(0,y,z) == bg ) maux(0,y,z) = 0;
        else                     maux(0,y,z) = INFTY;
      
        // Varredura para frente
        for( int x=1; x<nx; x++){
            if( IMG(x,y,z) == bg ) maux(x,y,z) = 0;
            else                     maux(x,y,z) = sum( 1, maux(x-1,y,z));   
        }
      
        // Varredura para trás
        for(int x=nx-2; x>=0; x--){    
            if( maux(x+1,y,z) < maux(x,y,z) ) 
                maux(x,y,z) = sum(1, maux(x+1,y,z));
        }
    }}

    // ---------------------------------------------------------------------------
    // PASSO 2: Varredura ao longo do eixo Y
    MTvi s(ny > nz ? ny : nz); // Reutiliza o vetor s para o maior entre ny e nz
    MTvi t(ny > nz ? ny : nz); // Reutiliza o vetor t
    int q, w;

    for( int x=0; x<nx; x++){
    for( int z=0; z<nz; z++){
        q=0;
        s[0] = 0;
        t[0] = 0;
    
        // Varredura para frente
        for( int u=1; u<ny; u++){ // u aqui representa o índice y
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
    
        // Varredura para trás
        for( int u=ny-1; u>=0; --u ){
            maux2(x,u,z) = F( u, s[q], prod( maux(x,s[q],z), maux(x,s[q],z)) );    
            if( u==t[q] ) q--;
        }
    }}

    // ---------------------------------------------------------------------------
    // PASSO 3: Varredura ao longo do eixo Z
    for( int x=0; x<nx; x++){
    for( int y=0; y<ny; y++){
        q=0;
        s[0] = 0;
        t[0] = 0;
    
        // Varredura para frente
        for( int u=1; u<nz; u++){ // u aqui representa o índice z
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
      
        // Varredura para trás
        for( int u=nz-1; u>=0; --u){
            EDT(x,y,u) = F( u, s[q], maux2(x,y,s[q]) );       
            if( u==t[q] ) 
                q--;
        }
    }}
}







#endif // EUCLIDIAN_DISTANCE_TRANSFORM





















