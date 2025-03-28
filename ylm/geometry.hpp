//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//  Funções que aplicam transformações na geometria das rochas.
//________________________________________________________
//A.Z. - 12/13 => Criação
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef FMM_GEOMETRY_HPP
#define FMM_GEOMETRY_HPP


// Padrão do C++
#include <iostream>
using namespace std;


// Minhas bibliotecas
#include "meusTipos.hpp"


// Para matrizes 3D com Boost
#include "boost/multi_array.hpp"
using namespace boost;
typedef multi_array<int , 3> matrix;




//------------------------------------------------------------------------------
// DESCRICAO:
//   Coloca paredes em uma geometria, conforme a direção de invasão
//
// RECEBE:
//   m => Matriz onde colocar as paredes
//   S => Número para colocar em um pixel sólido
//   a => Eixo
//
// ATENCAO:
//   m é alterada!
//
// RETORNA:
//   m => Matriz é alterada dentro da função, então retorna por referência
void walls( matrix &m, MTci &S, MTcs &a ){
  MTci nx = m.shape()[0];
  MTci ny = m.shape()[1];
  MTci nz = m.shape()[2];
  
  if( a=="x" ){
    for( int x=0; x<nx; x++ ){
      
      if( nz>1 ){  // Esses ifs são necessários para o caso de só haver um
                   // plano nessa direção (caso 2D), não colocar paredes nesse
                   // plano, senão eu tamparia a imagem toda.
        for( int y=0; y<ny; y++ ){
          m[x][y][   0] = S;
          m[x][y][nz-1] = S;
        }
      }

      if( ny>1 ){
        for( int z=0; z<nz; z++ ){
          m[x][   0][z] = S;
          m[x][ny-1][z] = S;
        }
      }
    }
    
    
  }else if( a=="y" ){
    for( int y=0; y<ny; y++ ){
      
      if( nz>1 ){
        for( int x=0; x<nx; x++ ){
          m[x][y][   0] = S;
          m[x][y][nz-1] = S;
        }
      }
  
      if( nx>1 ){
        for( int z=0; z<nz; z++ ){
          m[   0][y][z] = S;
          m[nx-1][y][z] = S;
        }
      }
    }

  }else if( a=="z" ){
    for( int z=0; z<nz; z++ ){
      
      if( ny>1 ){
        for( int x=0; x<nx; x++ ){
          m[x][   0][z] = S;
          m[x][ny-1][z] = S;
        }
      }
  
      if( nx>1 ){
        for( int y=0; y<ny; y++ ){
          m[   0][y][z] = S;
          m[nx-1][y][z] = S;
        }
      }
    }

  }
    
}








//------------------------------------------------------------------------------
// DESCRICAO:
//   Coloca uma membrana semipermeável (1px) em uma geometria em um valor
//   específico de y, geralmente 0 ou ny-1.
//
// RECEBE:
//   m => Matriz onde colocar a membrana
//   S => Número para colocar em um pixel sólido
//   P => Número para colocar em um pixel poro
//   p => Posição onde colocar a membrana
//   a => Eixo
//
// ATENCAO:
//   m é alterada!
//
// RETORNA:
//   m => Matriz é alterada dentro da função, então retorna por referência
void membrane( matrix &m, MTci &S, MTci &fill, MTci &p, MTcs &a ){
  MTci nx = m.shape()[0];
  MTci ny = m.shape()[1];
  MTci nz = m.shape()[2];

  if      ( a=="x" ){
    for( int z=0; z<nz; z++  ){
      
      MTci y0s = z%2;           // Início dos pontos sólidos, 0 ou 1
      MTci y0p = (y0s==0)? 1:0; // Início dos pontos porosos, 0 ou 1
      
      for( int y=y0s; y<ny; y+=2 ) m[p][y][z] = S;
      for( int y=y0p; y<ny; y+=2 ) m[p][y][z] = fill;
    }
        
  }else if( a=="y" ){
    for( int z=0; z<nz; z++  ){
      
      MTci x0s = z%2;           // Início dos pontos sólidos, 0 ou 1
      MTci x0p = (x0s==0)? 1:0; // Início dos pontos porosos, 0 ou 1
      
      for( int x=x0s; x<nx; x+=2 ) m[x][p][z] = S;
      for( int x=x0p; x<nx; x+=2 ) m[x][p][z] = fill;
    }
  }else if( a=="z" ){
    for( int y=0; y<ny; y++  ){
      
      MTci x0s = y%2;           // Início dos pontos sólidos, 0 ou 1
      MTci x0p = (x0s==0)? 1:0; // Início dos pontos porosos, 0 ou 1
      
      for( int x=x0s; x<nx; x+=2 ) m[x][y][p] = S;
      for( int x=x0p; x<nx; x+=2 ) m[x][y][p] = fill;
    }
  }

}



void surround( matrix &m, MTci &I){
  MTci nx = m.shape()[0];
  MTci ny = m.shape()[1];
  MTci nz = m.shape()[2];
  

    for( int x=0; x<nx; x++ ){
      
      if( nz>1 ){  // Esses ifs são necessários para o caso de só haver um
                   // plano nessa direção (caso 2D), não colocar paredes nesse
                   // plano, senão eu tamparia a imagem toda.
        for( int y=0; y<ny; y++ ){
          m[x][y][   0] = I;
          m[x][y][nz-1] = I;
        }
      }

      if( ny>1 ){
        for( int z=0; z<nz; z++ ){
          m[x][   0][z] = I;
          m[x][ny-1][z] = I;
        }
      }
    }
    
    

    for( int y=0; y<ny; y++ ){
      
      if( nz>1 ){
        for( int x=0; x<nx; x++ ){
          m[x][y][   0] = I;
          m[x][y][nz-1] = I;
        }
      }
  
      if( nx>1 ){
        for( int z=0; z<nz; z++ ){
          m[   0][y][z] = I;
          m[nx-1][y][z] = I;
        }
      }
    }


    for( int z=0; z<nz; z++ ){
      
      if( ny>1 ){
        for( int x=0; x<nx; x++ ){
          m[x][   0][z] = I;
          m[x][ny-1][z] = I;
        }
      }
  
      if( nx>1 ){
        for( int y=0; y<ny; y++ ){
          m[   0][y][z] = I;
          m[nx-1][y][z] = I;
        }
      }
    

  }
    
}



#endif // FMM_GEOMETRY_HPP





















