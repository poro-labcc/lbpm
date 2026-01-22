// //=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// // Classe Young_Laplace_Method, 2D e 3D
// //   
// //   Aplica o Young Laplace Method usado operações de Morfologia Matemática.
// //   Baseado no paper de Magnani et al, 2000 e bibliografias diversas de 
// //  Morfologia Matemática.
// //
// //
// //  ATENCAO ATENCAO ATENCAO
// //   Em todas as chamadas de função, para localizar o pixel uso a notação
// //   (x,y,z). Portanto, (coluna,linha,plano).
// //
// //=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// #ifndef YOUNG_LAPLACE_METHOD_HPP
// #define YOUNG_LAPLACE_METHOD_HPP


// // Padrão do C++
// #include <fstream>
// #include <iomanip>
// #include <stdint.h>
// using namespace std;



// // Minhas bibliotecas para o YLM
// #include "edt.hpp"
// #include "geometry.hpp"
// #include "component_labeling.hpp"


// // Minhas bibliotecas de uso geral
// #include "file_uti.hpp"
// #include "true_time.hpp"
// #include "meusTipos.hpp"
// #include "dataFileReader.hpp"

// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include "../common/Array.h"
// #include "../common/Domain.h"
// #include "../analysis/distance.h"
// // #include "../analysis/morphology.h"

// // Para matrizes 2D e 3D com Boost
// #include "boost/multi_array.hpp"
// typedef multi_array<int , 3> matrix;
// typedef multi_array<bool, 3> matrix_bool;

// // Desabilito as checagens do Boost para ganhar velocidade
// // Mas se houver algum problema, só vai dar pau no programa, nenhum aviso ;-)
// #define BOOST_DISABLE_ASSERTS




// // Matriz 2D STL para armazenar equivalências
// typedef vector< vector<int> > MTvvi;






// // -----------------------------------------------------------------------------
// // Classe
// class Young_Laplace_Method{
//   public:
  
//     Young_Laplace_Method( int, char *[] );
    
//     int Ndiam (void) const { return _d.size(); }
//     int diameter( MTci & );
  
//     void calc( MTci & );
    
    
//   public:
  
//     //int  _nthreads;             // Número de threads para usar no OpenMP
  
//     int _ny, _nx, _nz;          // Dimensões do Micromodelo
//     int dimy, dimx, dimz;       // salvo apenas os parâmetros de entrada; a ser utilizado no print .dat ADICIONADO AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
//     matrix _mmorig;             // Micromodelo original, não editável
//     matrix _mm, _mmcopy;                 // Micromodelo de trabalho, editável
    
//     int _NP;                    // Número de pixeis na parte porosa
    
//     MTvi _d;                    // Diâmetros (px). Só podem ser números ímpares
  
//     bool _wet, _compressible;   // Molhabilidade e compressibilidade
    
    
//     bool _iX, _iY, _iZ;         // Direção da invasão: x, y ou z
//     bool _iP, _iM;              // Sentido da invasão: + ou -
//     int _chamber_i, _chamber_o; // Onde estão os reservatórios?
    
    
//     int _memb;                  // Coloca membrana?
//     bool _wall;                 // Paredes?  
    
  
//     // Partes dos nomes das imagens de saída e do arquivo de dados e de tempo
//     string _outImgRoot, _whichImg, _outImgDir; 
//     ofstream _dat, _ftime;
    
    
//     True_Time _ttime;            // Medição do tempo
    
  
//     // Cores para as várias regiões da imagem.
//     int _I,   // I -> Fluido que está entrado
//         _O,   // O -> Fluido que está saindo
//         _S,   // S -> Parte sólida
//         _P;   // P -> Poros: estão preenchidos com O e serão invadidos por I
    
  
//     // Matriz booleanas para dizer se os pixeis são
//     matrix_bool _trapped; // regiões presas na invasão ou não
  
  
//     // Às vezes o fluido invasor não está só no reservatório, mas também
//     // fora, na imagem. Nesse caso, é preciso um loop a mais no cálculo
//     bool has_out_inlet;
  
  
//     // Matrizes auxiliares do tamanho da imagem.
//     // Não guardam nenhuma informação permanente
//     matrix _matrix1, _matrix2; // Matrizes usadas dentro da função calc
  
//     // Matrizes usadas EXCLUSIVAMENTE pela função euclidian_distance_transform
//     // São passadas por referência para evitar o overhead de serem criadas 
//     // e inicializadas toda hora.  
//     matrix _mx_edt1, _mx_edt2; 
    
    
//     // Vetores usados EXCLUSIVAMENTE pela função Component Labeling
//     // São passados por referência para evitar o overhead de serem criados 
//     // e inicializados toda hora.    
//     MTvi _next, _tail, _rtable;
  
    
//     // Transformada de Distância Euclidiana (EDT)
//     matrix _edt;              // EDT da imagem original
//   };
  
  
  
  
  
  
  
  
  
  
//   //------------------------------------------------------------------------------
//   // DESCRICAO:
//   //   Construtor
//   //   Lê arquivo de configurações e arquivo de imagem.
//   //   Inicializa parâmetros da classe para a simulação.
//   // RECEBE:
//   //   config => Nome do arquivo de configurações
//   Young_Laplace_Method::Young_Laplace_Method( int argc, char *argv[] ){
//     // string saux;


//     char LocalRankFilename[40];
//     string filename;
//     double Rcrit_new;
//     filename=argv[1];
//     Rcrit_new=0.f; 
//     NULL_USE( Rcrit_new );
//     // read the input database 
// 		auto db = std::make_shared<Database>( filename );
// 		auto domain_db = db->getDatabase( "Domain" );
//     auto ylm_db = db->getDatabase( "YLM" );

//     auto size = domain_db->getVector<int>( "N" );
//     _nx = size[0];
//     _ny = size[1];
//     _nz = size[2];

//     auto VoxelLabels = domain_db->getVector<int>( "VoxelLabels" );
//     _S = VoxelLabels[0];
//     _O = VoxelLabels[1];
//     _I = VoxelLabels[2];
//     _P = VoxelLabels[3];

//     auto READFILE = domain_db->getScalar<std::string>( "Filename" );
//     MTcs mmfile(READFILE);

//     _outImgRoot = ylm_db->getScalar<std::string>( "ImageRoot" );
//     _outImgDir = ylm_db->getScalar<std::string>( "ImageDir" );
//     _whichImg = ylm_db->getScalar<std::string>( "WhichImage" );
//     string memb = "none";
//     MTcs has_membrane(memb);
//     _wall = false;
//     _wet = ylm_db->getScalar<bool>( "Wetting" );
//     _compressible = false;


//     auto direction = ylm_db->getVector<int>( "Direction" );
//     _iX=false; _iY=false; _iZ=false;
//     _iM = false; _iP = false;
//     if(direction[0] != 0 && direction[1] == 0 && direction[2] ==0 ){_iX=true;if(direction[0] > 0){_iP = true;}else {_iM = true;}}
//     else if(direction[1] != 0 && direction[0] == 0 && direction[2] ==0 ){_iY=true;if(direction[1] > 0){_iP = true;}else {_iM = true;}}
//     else if(direction[2] != 0 && direction[0] == 0 && direction[1] ==0 ){_iZ=true;if(direction[2] > 0){_iP = true;}else {_iM = true;}}
//     else {aborta("Unkonwn direction.");}
    
//     auto diameters = ylm_db->getVector<int>( "Diameters" );
//     for (int dd = diameters[0]; dd <= diameters[1]; dd+= diameters[2]){
//       _d.push_back(dd);
//     }
//        // Ordena
//        sort( _d.begin(), _d.end() );
    
//        // Tira repetidos
//        MTvi::iterator it = unique( _d.begin(), _d.end() );
//        _d.resize( distance( _d.begin(),it ) );
     
//        // Aborta se não sobrou nada
//        if( _d.size()==0 )
//        aborta("It was impossible to create diameters array.");
     
//        // Se for não-molhante, reverte ordem do vetor
//        if( !_wet )  reverse( _d.begin(), _d.end() );
    
  
//     // Cria diretório para guardar imagens e arquivo de dados e de tempo
//     _outImgDir += "/";
//     mymkdir( _outImgDir );
  
    
   
//     // // ---------------------------------------------------------------------------
//     // // Lê a imagem
//     // //   Coloco duas linhas a mais. Uma em cima para representar o reservatório de
//     // // entrada de fluido e outra embaixo para o reservatório de fluido expulso.
//     // //   Todo o algoritmo segue normalmente com a imagem aumentada. Na hora de
//     // // gravar a imagem, retiro as duas linhas extras.
//     // // ---------------------------------------------------------------------------
//     double tt_read = _ttime(); // Tempo de leitura dos dados
  
  
//     // A leitura do arquivo exige pular os planos extra colocados
//     int x0=0, y0=0, z0=0;
//     int xM=_nx, yM=_ny, zM=_nz;
//     dimx = _nx; dimy = _ny; dimz = _nz; //ADICIONADO AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  
  

//     if(!_compressible){
//     // ---------------------------------------------------------------------------
//     // Inicio as variáveis para os reservatórios
//     // Coloco dois planos a mais, para conter os dois reservatórios
//     if     ( _iX ){ _nx += 2;   x0++;   xM = _nx-1; }
//     else if( _iY ){ _ny += 2;   y0++;   yM = _ny-1; }
//     else if( _iZ ){ _nz += 2;   z0++;   zM = _nz-1; }
//     // ---------------------------------------------------------------------------
//     }
  
//     // ---------------------------------------------------------------------------
//     // Membrana
//     // Quando tem membrana, coloco um plano a mais, entre o reservatório e a rocha
//     _memb = -1;
//     string pos_memb = has_membrane;
//     if( _iM ){
//       if     ( pos_memb=="begin" ) pos_memb="end";
//       else if( pos_memb=="end"   ) pos_memb="begin";
//     }
      
    
//     if( _iX ){
//       if( pos_memb == "begin" ){
//         _nx++;
//         _memb=1;
//         x0++; xM++;
//       }else if( pos_memb == "end" ){
//         _nx++;
//         _memb=_nx-2;
//         xM=_memb;
//       }
      
//     }else if( _iY ){
//       if( pos_memb == "begin" ){
//         _ny++;
//         _memb=1;
//         y0++; yM++;
//       }else if( pos_memb == "end" ){
//         _ny++;
//         _memb=_ny-2;
//         yM=_memb;
//       }
    
//     }else if( _iZ ){
//       if( pos_memb == "begin" ){
//         _nz++;
//         _memb=1;
//         z0++; zM++;
//       }else if( pos_memb == "end" ){
//         _nz++;
//         _memb=_nz-2;
//         zM=_memb;
//       }
    
//     }
//     // ---------------------------------------------------------------------------
  
  
  
  
//     // ---------------------------------------------------------------------------
//     // Se a imagem tem paredes, preciso reservar um espaço para elas também
//     // Se for uma imagem 2D, não coloco paredes no plano que só tem 1 px de largura
//     //
//     if( _wall ){
//       if( _iX ){
//         if( _ny>1 ){ _ny+=2;  y0++;  yM++; }
//         if( _nz>1 ){ _nz+=2;  z0++;  zM++; }
      
//       }else if( _iY ){
//         if( _nx>1 ){ _nx+=2;  x0++;  xM++; }
//         if( _nz>1 ){ _nz+=2;  z0++;  zM++; }
      
//       }else if( _iZ ){
//         if( _nx>1 ){ _nx+=2;  x0++;  xM++; }
//         if( _ny>1 ){ _ny+=2;  y0++;  yM++; }
        
//       }
//     }
//     // ---------------------------------------------------------------------------
  

//     if(_compressible) {
//       if( _nx>1 ){ _nx+=2;  x0++;  xM++; }
//       if( _ny>1 ){ _ny+=2;  y0++;  yM++; }
//       if( _nz>1 ){ _nz+=2;  z0++;  zM++; }
//     }
  
  
//     // ---------------------------------------------------------------------------
//     // A posição dos revervatórios depende do sentido de invasão
//     int iaux=-1;
//     if     ( _iX ){ iaux=_nx-1; }
//     else if( _iY ){ iaux=_ny-1; }
//     else if( _iZ ){ iaux=_nz-1; }
  
//     if( _iP ){
//       _chamber_i =    0;
//       _chamber_o = iaux;
//     }else if( _iM ){
//       _chamber_i = iaux;
//       _chamber_o =    0;
//     }
//     // ---------------------------------------------------------------------------
   
  
  
  
//     // ---------------------------------------------------------------------------
//     // As matrizes já foram inicializadas
//     // Preciso redimensionar para não dar problema
//     matrix::extent_gen extents;
//     _mmorig.resize  ( extents[_nx][_ny][_nz] );
//     _mm.resize      ( extents[_nx][_ny][_nz] );
//     _mmcopy.resize  ( extents[_nx][_ny][_nz] );
//     _edt.resize     ( extents[_nx][_ny][_nz] );
//     _matrix1.resize ( extents[_nx][_ny][_nz] );
//     _matrix2.resize ( extents[_nx][_ny][_nz] );
//     _mx_edt1.resize ( extents[_nx][_ny][_nz] );
//     _mx_edt2.resize ( extents[_nx][_ny][_nz] );
  
//     matrix_bool::extent_gen extentsb;
//     _trapped.resize( extentsb[_nx][_ny][_nz] );
  
  
//     // Quando a conexão é só entre primeiros vizinhos, podemos ter metade dos
//     // pixeis sendo desconectados. Pense num tabuleiro de xadrez 3D.
//     MTci N3 = _nx*_ny*_nz;
//     MTci max_eq = static_cast<int> ( ceil( static_cast<double>(N3)/2.0 ) );
//     _tail.resize( max_eq );
//     _next.resize( max_eq );
//     _rtable.resize( max_eq );
//     // ---------------------------------------------------------------------------
   
//   if(_compressible) surround(_mm, _I);
//   cout << _nx << _ny << _nz << endl;
  
  
//     // ---------------------------------------------------------------------------
//     // Coloco a membrana e depois a parede
//     if( _memb >0 ){
//       //MTci fill = (_wet)? _I:_P;
//       MTci fill = _P;
      
//       if     ( _iX ){ membrane( _mm, _S, fill, _memb, "x" ); }
//       else if( _iY ){ membrane( _mm, _S, fill, _memb, "y" ); }
//       else if( _iZ ){ membrane( _mm, _S, fill, _memb, "z" ); }    
//     }

//     if( _wall ){
//       if     ( _iX ){ walls( _mm, _S, "x" ); }
//       else if( _iY ){ walls( _mm, _S, "y" ); }
//       else if( _iZ ){ walls( _mm, _S, "z" ); }
//     }
//     // ---------------------------------------------------------------------------
  
    
  
  
  
  
//     // ---------------------------------------------------------------------------
//     // Inicializa o Micromodelo
    
//     // Coloca dois planos com Reservatórios
//     if(!_compressible){
//     // Invasão em x
//     if( _iX ){
//       for( int z=0; z<_nz; z++ ){
//       for( int y=0; y<_ny; y++ ){
//         _trapped[_chamber_i][y][z] = false;
//         _trapped[_chamber_o][y][z] = false;
        
//         _mm[_chamber_i][y][z] =_I;
//         _mm[_chamber_o][y][z] =_O;
//       }}
      
//     // Invasão em y    
//     }else if( _iY ){
//       for( int z=0; z<_nz; z++ ){
//       for( int x=0; x<_nx; x++ ){
//         _trapped[x][_chamber_i][z] = false;
//         _trapped[x][_chamber_o][z] = false;
        
//         _mm[x][_chamber_i][z] =_I;
//         _mm[x][_chamber_o][z] =_O;
//       }}
      
//     // Invasão em z
//     }else if( _iZ ){
//       for( int y=0; y<_ny; y++ ){
//       for( int x=0; x<_nx; x++ ){
//         _trapped[x][y][_chamber_i] = false;
//         _trapped[x][y][_chamber_o] = false;
        
//         _mm[x][y][_chamber_i] = _I;
//         _mm[x][y][_chamber_o] = _O;
//       }}
      
//     } }
      

    
  
//     //cout << "_nz = " << _nz << "; " << "z = " << z0 << " " << zM << endl;
//     //cout << "_ny = " << _ny << "; " << "y = " << y0 << " " << yM << endl;
//     //cout << "_nx = " << _nx << "; " << "x = " << x0 << " " << xM << endl;
  
  
  
//     // Lê em formato raw e já escreve em vtk
//     ifstream FRAW( mmfile.c_str() );
//     abriu( FRAW, mmfile );
//     unsigned char auxraw;
//     int pc;
   
  
//     // Imagem realmente lida do arquivo, sem os reservatórios
//     has_out_inlet=false;
//     for( int z=z0; z<zM; z++ ){
//     for( int y=y0; y<yM; y++ ){
//     for( int x=x0; x<xM; x++ ){
  
//       _trapped[x][y][z] = false;
//       FRAW >> auxraw;
//       pc = static_cast<int>( auxraw );
      
//       if( pc == _I )  has_out_inlet=true;
//       if( pc != _S )  _NP++;
  
//       if( pc!=_I  &&  pc!=_S  &&  pc!=_P  &&  pc!=_O ){  
//         aborta("Cor desconhecida em (" + ntos(x) + ", " +ntos(y) + ", " +ntos(z) + ")." );
//       }
      
//       _mm[x][y][z]     = pc;
//     }
//     }}
  
//     // Copia _mm em _mmorig
//     // Conta quantos pixeis porosos há, mas conta os reservatórios junto ...
//     _NP=0;
//     for( int z=0; z<_nz; z++ ){
//     for( int y=0; y<_ny; y++ ){
//     for( int x=0; x<_nx; x++ ){
//       iaux = _mm[x][y][z];
//       if( iaux!=_S ) _NP++;
//       _mmorig[x][y][z] = iaux;
//     }}}  
  
//     if(!_compressible){
//     // Desconta o tamanho dos reservatórios ...
//     if     ( _iX ){ _NP -= 2*_ny*_nz;  }
//     else if( _iY ){ _NP -= 2*_nx*_nz;  }
//     else if( _iZ ){ _NP -= 2*_nx*_ny;  } }
//     else{
//       _NP = _NP - 2*(_nx*_ny+_nx*_nz+_ny*_nz) + (4*(_nx+_ny+_nz) - 8);
//     }
  
  
  
    
//     tt_read = _ttime() - tt_read; // Tempo de leitura dos dados
//     // ---------------------------------------------------------------------------
  
   
    
//     // ---------------------------------------------------------------------------
//     // Calcula Transformada de Distância Euclidiana da imagem _mmorig e coloca o
//     // resultado na matriz _edt. Usa _S como cor de fundo
//     double tt_edt = _ttime(); // Tempo para calcular a edt
//     // euclidian_distance_transform( _S, _mmorig, _edt, _nthreads, _mx_edt1, _mx_edt2 );
//     euclidian_distance_transform( _S, _mmorig, _edt, _mx_edt1, _mx_edt2 );
//     tt_edt = _ttime() - tt_edt; // Tempo para calcular a edt
//     // ---------------------------------------------------------------------------
  
  
  
//     // ---------------------------------------------------------------------------
//     // Arquivo de saída com resultados da invasão
//     MTvs vec(6), cmt(4);
    
//     vec[0] = "Step";
//     vec[1] = "Diameter (px)";
//     vec[2] = "Number of pixels occupied by inlet fluid.";
//     vec[3] = "Number of pixels occupied by inlet fluid / Number of porous pixels";
//     vec[4] = "Number of pixels occupied by outlet fluid.";
//     vec[5] = "Number of pixels occupied by outlet fluid  / Number of porous pixels";
    
//     if(true) {
  
//       cmt[0] = "Image width  (x) (px)     : " + ntos( dimx );
//       cmt[1] = "Image height (y) (px)     : " + ntos( dimy );
//       cmt[2] = "Image planes (z) (px)     : " + ntos( dimz );
//       cmt[3] = "Number of porous pixels   : " + ntos( _NP  ); }
  
//     // else {
  
//     //   int camadas = 2*(dimx*dimy+dimx*dimz+dimy*dimz) - (4*(dimx+dimy+dimz) - 8);
//     //   int N_NP = _NP - camadas;
  
//     //   cmt[0] = "Image width  (x) (px)     : " + ntos( dimx - 2) + "       (disconsidering surrounding layers)";
//     //   cmt[1] = "Image height (y) (px)     : " + ntos( dimy - 2) + "       (disconsidering surrounding layers)";
//     //   cmt[2] = "Image planes (z) (px)     : " + ntos( dimz - 2) + "       (disconsidering surrounding layers)";
//     //   cmt[3] = "Number of porous pixels   : " + ntos(   N_NP  ) + "       (disconsidering surrounding layers)";
  
//     // }
  
//     // Cria arquivo de saída
//     // Deixo que o destrutor implícito feche o arquivo
//     string saux = _outImgDir + "/" + _outImgRoot + ".dat";
//     _dat.open( saux.c_str() );
//     abriu( _dat, saux );
//     outputFileHead( argc, argv, _dat, vec, cmt );  
//     // ---------------------------------------------------------------------------
  
  
  
//     // ---------------------------------------------------------------------------
//     // Arquivo de tempos
//     saux = _outImgDir + "/" + _outImgRoot + "_time.txt";
//     _ftime.open( saux.c_str() );
//     _ftime << scientific << setprecision(6);
    
//     abriu( _ftime, saux );
//     _ftime << "# Arquivo com tempos de execução do código" << endl;
//     _ftime << "# Tempos em segundos" << endl;
//     _ftime << endl;
    
//     _ftime << "# " << tt_read << " # Tempo para ler e inicializar os dados" << endl;
//     _ftime << "# " << tt_edt  << " # Tempo para calcular a EDT na inicialização" << endl;
    
    
//     _ftime << "\n\n"
//            << "# A partir daqui, listo os tempos que cada passo levar para executar vários pedaços do programa.\n"
//            << "# \n"
//            << "# Col  1: Número do passo\n"
//            << "# Col  2: Diâmetro do passo\n"
//            << "# Col  3: tt_step\n"
//            << "# Col  4: tt_loop1\n"
//            << "# Col  5: tt_edt1\n"
//            << "# Col  6: tt_loop2\n"
//            << "# Col  7: tt_cl1\n"
//            << "# Col  8: tt_loop3\n"
//            << "# Col  9: tt_loop1w\n"
//            << "# Col 10: tt_cl1w\n"
//            << "# Col 11: tt_loop2w\n"
//            << "# Col 12: tt_loop4\n"
//            << "# Col 13: tt_loop5\n"
//            << "# Col 14: tt_loop1c\n"
//            << "# Col 15: tt_cl1c\n"
//            << "# Col 16: tt_loop2c\n"
//            << "# Col 17: tt_end\n"
//            << endl;
//     // ---------------------------------------------------------------------------  
//   }
  
  
  
  
  
  
  
  
  
  
  
//   //------------------------------------------------------------------------------
//   // DESCRICAO:
//   //   Retorna o valor de um diâmetro para um dado passo
//   // RECEBE:
//   //   step => Número do passo atual
//   int Young_Laplace_Method::diameter( MTci &step ){
//     if( step<0 || step>=_d.size() )
//       aborta( "Valor inválido para o passo." );
//     return _d[step];
//   }
  
  
  
  
  
  
  
  
  
//   //------------------------------------------------------------------------------
//   // DESCRICAO:
//   //   Calcula a invasão para um diâmetro específico
//   //
//   //   Por uma questão de eficiência, não as etapas de modo sequencial, mas
//   //   aproveito os loops sobre a imagem para fazer pedaços de várias etapas.
//   //   Isso torna o código bem mais confuso, mas economiza loops, o que torna
//   //   o programa mais eficiente. Não sei se o custo da confusão paga o ganho
//   //   em eficiência.
//   // 
//   // RECEBE:
//   //   step => Número do passo para calcular
//   void Young_Laplace_Method::calc( MTci &step ){
//     int iaux;
    
    
//     MTci D = this->diameter(step);
//     _ftime << setw(3) << step << setw(4) << D;
  
  
//     double tt_step = _ttime(); // Tempo para executar esse passo inteiro
  
    
    
//     //----------------------------------------------------------------------------
//     // Na erosão, preciso eliminar os pixeis que não fazem parte da imagem erodida
//     // Seja (w,h) o centro de uma bola e (x,y) um pixel sólido mais próximo.
//     // Temos:
//     //   1 - EDT = (w-x)^2 + (h-y)^2
//     //   2 - O pixel (x,y) faz parte da bola de diâmetro D e centro (w,h) se
//     //       (w-x)^2 + (h-y)^2 <= D^2/4
//     //   3 - Na erosão, quero tirar todos os pixeis que não fazem parte da bola,
//     //       então, excluo os pixeis com EDT^2 que não satisfazem o critério 2.
//     //
//     //   Logo, comparo EDT com D^2/4. Simples assim.
//     //
//     MTcd D24 = D*D/4.0;
//     //----------------------------------------------------------------------------
  
  
//     // Índices para Background e Foreground
//     // É necessário que tenham ESTES valores
//     MTci B=0, F=1;
  
  
//     //----------------------------------------------------------------------------
//     // Realiza a operação de abertura para obter a região H do Magnani (eq. 5)
//     //
//     // A abertura é igual a uma erosão seguida de uma dilatação.
//     //
//     // Aplico o seguinte algoritmo baseado na Transformada de Distância Eucliana:
//     //   It is then straightforward to perform an erosion with a disc of radius 
//     //   r simply by removing all pixels whose distance label is less than r. A
//     //   dilation is similarly performed by eroding the background.
//     // descrito no livro IMAGE ANALYSIS FOR THE BIOLOGICAL SCIENCES, de
//     // C A GLASBEY and G W HORGAN, capítulo 5, página 10.
//     //----------------------------------------------------------------------------
  
  
//     //----------------------------------------------------------------------------
//     // LOOP 1:
//     //   * Determinar a região H por abertura:
//     //     - 1ª erosão -> matrix1
//     //   * Iguala _mm a _mmorig
//     //
//     double tt_loop1 = _ttime();
//     // #pragma omp parallel for num_threads (_nthreads)
//     for( int x=0; x<_nx; x++ ){
//     for( int y=0; y<_ny; y++ ){
//     for( int z=0; z<_nz; z++ ){
    
//       // Em _matrix1 coloco o resultado da 1ª erosão, que consiste em marcar
//       // os pixeis com r maior do que o do passo como Foreground e todo resto 
//       // como Background.
  
//       //   Além da região da erosão, preciso incluir os reservatórios para
//       // poder fazer o Component Labeling.
//       //   Depois eu tiro os pixeis do reservatório antes de fazer a dilatação,
//       // senão a minha operação não seria uma abertura verdadeira, pois eu
//       // dilataria regiões que não eram resultado da erosão, mas partes do
//       // reservatório.
//       //   Entretanto, não posso fazer o Component Labeling sem incluir os
//       // reservatórios, porque o reservatório pode ter várias partes desconectadas
//       // depois da erosão. Incluí-los conecta tudo, mas cria uma artefato que
//       // precisa ser corrigido antes da dilatação.
//       // 
//       //   Guardo a informação desses pixeis a mais colocados na erosão em
//       // matrix2, para poder tirá-los depois
//       _matrix2[x][y][z] = F;
      
//       // A região erodida é só aquela que EDT >= D24
//       if( _edt[x][y][z]>=D24 ){
//         _matrix1[x][y][z] = F;
        
      
//       }else{
//          _matrix1[x][y][z] = B;
        
  
//         // Aqui, tudo deveria ser background, mas tenho o reservatório 
//         // para incluir, por causa do que falei acima.
//         // Mas só está funcionando para líquidos não molhantes por enquanto!
//         if( !_wet && _mmorig[x][y][z]==_I ){
//            _matrix1[x][y][z] = F;
           
//            // Significará que foi colocado por ser do reservatório    
//            _matrix2[x][y][z] = B;  
//         }
//       }
      
      
//       // Se não apagar os resquícios da invasão anterior, dá erro no caso
//       // de ser molhante. Acho que é por que como já tem regiões invadidas,
//       // cria-se uma conexão com toda a imagem ... pensar melhor aqui.
//       // Mas quando não é molhante, se apago esses registros,
//       // acontece de regiões que já foram invadidas voltarem a ser não invadidas
//       // com um aumento do raio. Isso acontece por que as esferas não são
//       // realmente esféricas.
//       if( _wet )
//         _mm[x][y][z] = _mmorig[x][y][z];
      
//     }}}
//     tt_loop1 = _ttime() - tt_loop1;
//     //----------------------------------------------------------------------------
  
  
  
  
  
//     //----------------------------------------------------------------------------
//     // Encontro as regiões conexas
//     int xaux=_nx/2, yaux=_ny/2, zaux=_nz/2;
//     if     ( _iX ){ xaux=_chamber_i; }
//     else if( _iY ){ yaux=_chamber_i; }
//     else if( _iZ ){ zaux=_chamber_i; }
//     int chamber_label=0;
    
//     if( !_wet ){
//       // component_labeling( _matrix1, F, B, _next, _tail, _rtable, _nthreads );
//       component_labeling( _matrix1, F, B, _next, _tail, _rtable );
    
//       // Pego o label do reservatório  
//       chamber_label =  _matrix1[xaux][yaux][zaux];
      
//       // Tudo que não for conectado ao reservatório é colocado como background
//       // Os pixeis que eram só do reservatório, não da região erodida, também
//       // são marcados como background
//       // #pragma omp parallel for num_threads (_nthreads)
//       for( int x=0; x<_nx; x++ ){
//       for( int y=0; y<_ny; y++ ){
//       for( int z=0; z<_nz; z++ ){
//         if( _matrix1[x][y][z] == chamber_label ){
    
//           _matrix1[x][y][z] = _matrix2[x][y][z];
          
//         }else{
//           _matrix1[x][y][z] = B;
//         }
//       }}}  
//     }
//     //----------------------------------------------------------------------------
  
  
//     //for( int y=0; y<_ny; y++ ){
//     //for( int x=0; x<_nx; x++ ){
//       //cout << _matrix1[x][y][0] << " ";
//     //}cout << endl;}
  
  
  
  
  
//     // Agora, faço a EDT do negativo, dizendo que os pixeis Foreground são 
//     // Background 
//     // Calcula EDT da imagem _matrix1 e coloca o resultado em _matrix2.
//     double tt_edt1 = _ttime();
//     // euclidian_distance_transform( F, _matrix1, _matrix2, _nthreads, _mx_edt1, _mx_edt2  );
//     euclidian_distance_transform( F, _matrix1, _matrix2, _mx_edt1, _mx_edt2  );
//     tt_edt1 = _ttime() - tt_edt1;
  
  
  
//     //----------------------------------------------------------------------------
//     // LOOP 2a:
//     //   * Determinar a região H por abertura:
//     //     - 2ª erosão
//     //   * Marca as regiões ocupadas pelos fluidos na imagem _mm
//     double tt_loop2 = _ttime();
//     // #pragma omp parallel for num_threads (_nthreads)
//     for( int x=0; x<_nx; x++ ){
//     for( int y=0; y<_ny; y++ ){
//     for( int z=0; z<_nz; z++ ){
      
//       // Erosão do EDT negativo: se _matrix2 < D24, é background
//       // Portanto, é uma região ocupada por fluido invasor
//       // As outras regiões, se não for sólido, são ocupadas pelo fluido expulso
//       if( _matrix2[x][y][z] < D24 ){
//         _mm[x][y][z] = _I;
//       }else{
//        if( _mm[x][y][z] == _P )
//          _mm[x][y][z] = _O;          
//       }
     
//       // Prepara para encontrar componentes desconexas
//       _matrix1[x][y][z] = (_mm[x][y][z]==_I)? F:B;
//     }}}
//     tt_loop2 = _ttime() - tt_loop2;
    
    
    
//     //----------------------------------------------------------------------------
//     // LOOP 2b:
//     //   * Componentes desconexas do fluido invasor em matrix1
//     double tt_cl1 = _ttime();
//     // component_labeling( _matrix1, F, B, _next, _tail, _rtable, _nthreads );
//     component_labeling( _matrix1, F, B, _next, _tail, _rtable );
//     tt_cl1 = _ttime() - tt_cl1;
//     chamber_label =  _matrix1[xaux][yaux][zaux];
  
  
  
//     //----------------------------------------------------------------------------
//     // LOOP 3:
//     //   * Substitui índices equivalentes
//     //   * Determina região G
//     //   * Determina Omega = K( G, Bi )
//     //   * Determina que tipo de fluido haverá em cada pixel de acordo com Omega
//     //   * Colocar as cores da imagem final se for o momento oportuno
//     //
  
//     // !!!!!!!!!!!!!!!!!
//     // Não sei porque, mas não dá para paralelizar esse loop. Dá problema na
//     // embebição!!!!
//     //       
//     // TESTAR: É PORQUE iaux É SHARED E CADA THREAD MUDA O VALOR DELA INDEPENDENTEMENTE
//     //         SÓ PRECISO CRIAR UMA VARIÁVEL PRIVATE PARA ISSO E PRONTO:
//     //
//     // #pragma parallel for num_threads(_nthreads) private(iaux)
//     //
//     double tt_loop3 = _ttime();
//     for( int x=0; x<_nx; x++ ){
//     for( int y=0; y<_ny; y++ ){
//     for( int z=0; z<_nz; z++ ){
  
//       // Região G: L U H   eq.6
//       // Nesse passo, a região H é aquela em que tenho fluido invasor, pois
//       // depois da abertura coloquei fluido invasor em todo lugar onde cabia
//       // para fazer o Component Labeling
//       // A região L é somente a linha do reservatório do fluido invasor
//       bool rG=false;
//       if     ( _iX ){ rG  = (_mm[x][y][z]==_I || x==_chamber_i); }
//       else if( _iY ){ rG  = (_mm[x][y][z]==_I || y==_chamber_i); }
//       else if( _iZ ){ rG  = (_mm[x][y][z]==_I || z==_chamber_i); }
  
  
//       // Operador K, gero a região Omega   eq.11
//       // Não basta ser da região G, é preciso que esteja conectado ao reserv.
//       // Ou seja, precisa ter o mesmo label que os pixeis do reservatório
//       bool rO=false;
//       if( rG && _matrix1[x][y][z] == chamber_label )
//         rO = true;
  
//       // Coloco a cor certa, de acordo com o tipo de pixel
//       if( rO ){ 
//         _mm[x][y][z] = _I;
//       }else{
//         if( _mm[x][y][z] != _S )
//           _mm[x][y][z] = _O;
//       }
        
        
//       // Caso o fluido invasor seja molhante, preciso de alguns passos a mais,
//       // onde será necessário saber se um pixel é da região Omega e da região G
//       // Guardo essa informação em _matrix2      
//       if( _wet ){
//         iaux                = 0;  // pixel não é nem de G nem de O
//         if( rG ) iaux       = 1;  // pixel é só de G
//         if( rO ) iaux       = 2;  // pixel é só de O
//         if( rG && rO ) iaux = 3;  // pixel é de G e O
//         _matrix2[x][y][z] = iaux;
//       }
      
//     }}}
//     tt_loop3 = _ttime() - tt_loop3;
  
  
  
  
//     //----------------------------------------------------------------------------
//     // Caso seja um fluido molhante, são necessários alguns passos a mais:
//     //   1 - Calcular região X = Gc U Omega
//     //   2 - Calcular as componentes desconexas da região X
//     //   3 - Calcular K( X, Bi )
//     double tt_loop1w = 0;
//     double tt_cl1w = 0;
//     double tt_loop2w = 0;
  
//     if( _wet ){
  
//       //--------------------------------------------------------------------------
//       // LOOP w1a:
//       //   * Determina se um pixel é da região X = Gc U Omega1
//       tt_loop1w = _ttime();
//       // #pragma omp parallel for num_threads (_nthreads)
//       for( int x=0; x<_nx; x++ ){
//       for( int y=0; y<_ny; y++ ){
//       for( int z=0; z<_nz; z++ ){
  
//         // Pixel está em O1 e G?
//         iaux = _matrix2[x][y][z];
//         MTcb rG = iaux==1 || iaux==3;
//         MTcb rO = iaux==2 || iaux==3;
        
//         // Região Gc: F - G   eq.7
//         // Ou seja, o pixel está em F mas não está em G
//         MTcb rGc = _mm[x][y][z]!=_S &&  !rG;
        
//         // Nova região é composta por  Gc U O1
//         bool rX = false;
//         if( rGc || rO ) rX=true;
//         _matrix1[x][y][z] = (rX)? F:B;
//       }}}
//       tt_loop1w = _ttime() - tt_loop1w;
        
        
//       // Componentes desconexas da região X em matrix1
//       tt_cl1w = _ttime();
//       // component_labeling( _matrix1, F, B, _next, _tail, _rtable, _nthreads );
//       component_labeling( _matrix1, F, B, _next, _tail, _rtable);
//       tt_cl1w = _ttime() - tt_cl1w;
    
//       xaux=_nx/2, yaux=_ny/2, zaux=_nz/2;
//       if     ( _iX ){ xaux=_chamber_i; }
//       else if( _iY ){ yaux=_chamber_i; }
//       else if( _iZ ){ zaux=_chamber_i; }
//       chamber_label = _matrix1[xaux][yaux][zaux];
  
  
//       //--------------------------------------------------------------------------
//       // LOOP w2:
//       //   * Substitui índices equivalentes
//       //   * Determina Omega2 = K( X, Bi )
//       //   * Determina que tipo de fluido haverá em cada pixel de acordo com Omega2
//       //   * Colocar as cores da imagem final se for o momento oportuno
//       //
//       tt_loop2w = _ttime();
//       // #pragma omp parallel for num_threads (_nthreads)
//       for( int x=0; x<_nx; x++ ){
//       for( int y=0; y<_ny; y++ ){
//       for( int z=0; z<_nz; z++ ){
  
//         // Operador K, gero a região Omega2   eq.11
//         // Não basta ser da região X, é preciso que esteja conectado ao reserv.
//         // Ou seja, precisa ter o mesmo label que os pixeis do reservatório
//         bool rO2=false;
//         if( _matrix1[x][y][z] == chamber_label )
//           rO2 = true;
    
//         // Coloco a cor certa, de acordo com o tipo de pixel
//         if( rO2 ){ 
//           _mm[x][y][z] = _I;
//         }else{
//           if( _mm[x][y][z] != _S )
//             _mm[x][y][z] = _O;
//         }
//       }}}
//       tt_loop2w = _ttime() - tt_loop2w;
//     }
     
  
  
          
//     // Os reservatórios, sempre contém os fluidos originais, independente do
//     // processo de invasão atingí-los ou não.
//     double tt_loop4 = _ttime();
    
//     int x0=0  , y0=0  , z0=0;
//     int xM=_nx, yM=_ny, zM=_nz;
    
//     if( _iX ){
//       // #pragma omp parallel for num_threads (_nthreads)
//       for( int y=0; y<_ny; y++ ){
//       for( int z=0; z<_nz; z++ ){
//         _mm[_chamber_i][y][z] = _I;
//         _mm[_chamber_o][y][z] = _O;
//       }}
      
//       x0=1;
//       xM--;
      
//     // Invasão em y    
//     }else if( _iY ){
//       // #pragma omp parallel for num_threads (_nthreads)
//       for( int x=0; x<_nx; x++ ){
//       for( int z=0; z<_nz; z++ ){
//         _mm[x][_chamber_i][z] = _I;
//         _mm[x][_chamber_o][z] = _O;
//       }}
      
//       y0=1;
//       yM--;
      
//     // Invasão em z
//     }else if( _iZ ){
//       // #pragma omp parallel for num_threads (_nthreads)
//       for( int x=0; x<_nx; x++ ){
//       for( int y=0; y<_ny; y++ ){
//         _mm[x][y][_chamber_i] = _I;
//         _mm[x][y][_chamber_o] = _O;
//       }}
      
//       z0=1;
//       zM--;    
//     }
//     tt_loop4 = _ttime() - tt_loop4;      
     
  
     
//     // Se na imagem original o pixel era de um fluido invasor, agora ele pode
//     // ter sido substituido por um pixel expulso, pois poderia estar sem
//     // conexão com o reservatório. É o que acontece na experiência de intrusão.
//     // Fazemos uma drenagem, onde várias regiões ficam presas e depois uma
//     // embebição, onde essas regiões presas anteriormente não tem conexão com
//     // o novo reservatório invasor.
//     // Nesse caso, preciso colocar esses pixeis como o fluido invasor novamente
//     double tt_loop5 = _ttime();
//     if( has_out_inlet ){
//       // #pragma omp parallel for num_threads (_nthreads)
//       for( int x=x0; x<xM; x++ ){
//       for( int y=y0; y<yM; y++ ){
//       for( int z=z0; z<zM; z++ ){
//         if( _mmorig[x][y][z]==_I )
//           _mm[x][y][z] = _I;
//       }}}
//     }
//     tt_loop5 = _ttime() - tt_loop5;      
  
  
  
//     //----------------------------------------------------------------------------
//     // Caso seja um fluido incompressível, são necessários alguns passos a mais:
//     //   1 - Encontrar as regiões desconexas da imagem atual que não tem conexão
//     //       com o reservatório de fluido expulso
//     //   2 - Marcar os pixeis dessas regiões e mais os pixeis presos nos passos
//     //       anteriores como pixeis presos
//     //   3 - Colocar a cor desses pixeis como a cor do fluido expulso
//     double tt_loop1c=0;
//     double tt_cl1c=0;
//     double tt_loop2c=0;
    
//     if( !_compressible ){
  
//       //--------------------------------------------------------------------------
//       // LOOP c1a:
//       //   * Determina se um pixel é da região presa
//       tt_loop1c = _ttime();
//       // #pragma omp parallel for num_threads (_nthreads)
//       for( int x=0; x<_nx; x++ ){
//       for( int y=0; y<_ny; y++ ){
//       for( int z=0; z<_nz; z++ ){
  
//         // Se o pixel já estava preso no passo anterior, marco ele preso nesse
//         // também. É que parte das regiões presas no passo anterior foram
//         // invadidas nesse. Preciso "recuperar" esses pixeis como sendo presos.
//         if( _trapped[x][y][z] )
//           _mm[x][y][z] = _O;
  
//         _matrix1[x][y][z] = (_mm[x][y][z]==_O)? F:B;
//       }}}
//       tt_loop1c = _ttime() - tt_loop1c;      
      
      
      
//       //--------------------------------------------------------------------------
//       //   * Componentes desconexas da região presa em matrix1
//       tt_cl1c = _ttime();
//       // component_labeling( _matrix1, F, B, _next, _tail, _rtable, _nthreads );
//       component_labeling( _matrix1, F, B, _next, _tail, _rtable );
//       tt_cl1c = _ttime() - tt_cl1c;
      
    
//       xaux=_nx/2, yaux=_ny/2, zaux=_nz/2;
//       if     ( _iX ){ xaux=_chamber_o; }
//       else if( _iY ){ yaux=_chamber_o; }
//       else if( _iZ ){ zaux=_chamber_o; }
//       chamber_label = _matrix1[xaux][yaux][zaux];
      
      
//       //--------------------------------------------------------------------------
//       // LOOP c2:
//       //   * Substitui índices equivalentes
//       //   * Todos as regiões desconexas que não tem conexão com o reservatório
//       //     de fluido expulso são consideradas presas. Marco esses pixeis como
//       //     presos
      
//       // AQUI PRECISO MELHORAR
//       // É POSSÍVEL QUE O RESERVATÓRIO DE FLUIDO EXPULSO TENHA VÁRIAS REGIÕES
//       // DESCONEXAS? ENTÃO, O MAIS CORRETO SERIA FAZER UM LOOP PELO RESERVATÓRIO
//       // GUARDANDO OS LABELS DAS VÁRIAS REGIÕES. DEPOIS, QUALQUER REGIÃO
//       // QUE NÃO ESTEJA CONECTADA A UMA DESSAS REGIÕES SERÁ CONSIDERADA DESCONEXA
//       // E PORTANTO, PRESA
  
//       tt_loop2c = _ttime();
//       // #pragma omp parallel for num_threads (_nthreads)
//       for( int x=0; x<_nx; x++ ){
//       for( int y=0; y<_ny; y++ ){
//       for( int z=0; z<_nz; z++ ){    
  
//         // É uma região presa?
//         if( _mm[x][y][z]==_O &&  _matrix1[x][y][z] != chamber_label )
//           _trapped[x][y][z]=true;
//       }}}
//       tt_loop2c = _ttime() - tt_loop2c;
//     }
  
        
//     // ---------------------------------------------------------------------------
//     // Cria o arquivo que contém os dados da "imagem"
  
  
//     // Crio os arquivos de saída?
//     bool createRAW=false;
//     if( _whichImg=="all" ){
//       createRAW=true;
//     }else{
//       if( _whichImg=="last" && step==_d.size()-1 )
//         createRAW=true;
//     }
    
//     double tt_end = _ttime();
  
  
//     //AQUI SALVA OS ARQUIVOS (TIRAR RESERVATÓRIOS)


//     uint8_t aux_raw;
//     FILE *FRAW;
//     if( createRAW ){
      
//       xaux=_nx, yaux=_ny, zaux=_nz;
//       if     ( _iX ){ xaux = _nx-2; }
//       else if( _iY ){ yaux = _ny-2; }
//       else if( _iZ ){ zaux = _nz-2; }
      
//       int ndig = max(_d[0],_d[_d.size()-1]);
//       string saux = ".s" + ntos(ndig,step) + ".d" + ntos(ndig,D);
//       MTcs fraw = _outImgDir + "/" + _outImgRoot + saux + ".raw";
  
//       FRAW = fopen64(fraw.c_str(), "wb");
//     }
        
        
//     // ---------------------------------------------------------------------------
//     // Último loop para colocar as cores definitivas na imagem e fazer
//     // as contagens.
    
//     // Os loops não são completos em algumas dimensões, pois não contabilizo os
//     // reservatórios
//     x0=0  , y0=0  , z0=0;
//     xM=_nx, yM=_ny, zM=_nz;

//     if(!_compressible){
//     if     ( _iX ){ x0=1; xM = _nx-1; }
//     else if( _iY ){ y0=1; yM = _ny-1; }
//     else if( _iZ ){ z0=1; zM = _nz-1; } }
//     else{ x0=1; xM = _nx-1; y0=1; yM = _ny-1; z0=1; zM = _nz-1;}
    
  
//     // aux_raw = static_cast<uint8_t>(iaux);
//     // fwrite(&aux_raw, sizeof(uint8_t), 1, outfile);








  
    
//       //SALVAMENTO ANTIGO:

//     // int Ninlet=0, Noutlet=0, Nporo=0;
//     // for( int z=z0; z<zM; z++ ){
//     // for( int y=y0; y<yM; y++ ){
//     // for( int x=x0; x<xM; x++ ){    
      
//     //   iaux = _mm[x][y][z];
//     //   if     ( iaux==_I ) Ninlet++;
//     //   else if( iaux==_O ) Noutlet++;

//     //   if(iaux != _S ) Nporo++;
      
//     //   if( createRAW ){
//     //     aux_raw = static_cast<uint8_t>(iaux);
//     //     fwrite(&aux_raw, sizeof(uint8_t), 1, FRAW);
//     //   }
//     // }
//     // }
//     // }
//     // if( createRAW ) fclose(FRAW); 
//     // tt_end = _ttime() - tt_end;






//       //SALVEAMENTO NOVO (com aging):

//       int dx[26] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1, 0};
//       int dy[26] = {-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1};
//       int dz[26] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};


//     int Ninlet=0, Noutlet=0, Nporo=0;
//     for( int z=z0; z<zM; z++ ){
//     for( int y=y0; y<yM; y++ ){
//     for( int x=x0; x<xM; x++ ){    
      

//       iaux = _mm[x][y][z];
      
//       if (_mmcopy[x][y][z] != -_I && _mmcopy[x][y][z] != -_O)
//       _mmcopy[x][y][z] = iaux;
//       if (iaux == _I) {
//         Ninlet++;
//         // percorre 26 vizinhos
//         for(int n=0; n<26; n++){
//             int xn = x + dx[n];
//             int yn = y + dy[n];
//             int zn = z + dz[n];
//             // checa limites
//             if(xn>=x0 && xn<xM && yn>=y0 && yn<yM && zn>=z0 && zn<zM){
//                 if(_mm[xn][yn][zn] == _S) _mmcopy[xn][yn][zn] = -_I;
//               }
//           }
//       }
//       else if (iaux == _O) {
//         Noutlet++;
//         // percorre 26 vizinhos
//         for(int n=0; n<26; n++){
//             int xn = x + dx[n];
//             int yn = y + dy[n];
//             int zn = z + dz[n];
//             // checa limites
//             if(xn>=x0 && xn<xM && yn>=y0 && yn<yM && zn>=z0 && zn<zM){
//                 if(_mm[xn][yn][zn] == _S) _mmcopy[xn][yn][zn] = -_O;
//               }
//           }
//       }

//       if(iaux != _S ) Nporo++;
      
//     }
//     }
//     }

//     for( int z=z0; z<zM; z++ ){
//     for( int y=y0; y<yM; y++ ){
//     for( int x=x0; x<xM; x++ ){    
      
//       iaux = _mmcopy[x][y][z];

//       if( createRAW ){
//         aux_raw = static_cast<uint8_t>(iaux);
//         fwrite(&aux_raw, sizeof(uint8_t), 1, FRAW);
//       }}}}


//     if( createRAW ) fclose(FRAW); 
//     tt_end = _ttime() - tt_end;

    
 
  
//     //Modificação feita por: Thomas Carmo
//     // int camadas = 2*(dimx*dimy+dimx*dimz+dimy*dimz) - (4*(dimx+dimy+dimz) - 8); //ADICIONADO AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
//     // int N_NP = _NP - camadas;
  
//     // Salva contagens
  
//     // if(!_compressible){
//     if(true){
//     _dat << setprecision(6) << step << " " << D << " " << Ninlet << " " << Ninlet/(1.0*Nporo) << " " << Noutlet << " " << Noutlet/(1.0*Nporo) << endl;
//     }else{
   
//     // _dat << setprecision(6) << step << " " << D << " " << Ninlet-camadas << " " << (Ninlet-camadas)/((1.0*N_NP)) << " " << Noutlet << " " << Noutlet/((1.0*N_NP)) << endl;
//     }
//     // Salva tempos
//     tt_step = _ttime() - tt_step;  
//     _ftime << setw(14) << tt_step
//            << setw(14) << tt_loop1
//            << setw(14) << tt_edt1
//            << setw(14) << tt_loop2
//            << setw(14) << tt_cl1
//            << setw(14) << tt_loop3
//            << setw(14) << tt_loop1w
//            << setw(14) << tt_cl1w
//            << setw(14) << tt_loop2w
//            << setw(14) << tt_loop4
//            << setw(14) << tt_loop5
//            << setw(14) << tt_loop1c
//            << setw(14) << tt_cl1c
//            << setw(14) << tt_loop2c
//            << setw(14) << tt_end
//            << endl;
//   }
  
  
  
  
  
  
  
  
  
  
  
  
//   #endif // YOUNG_LAPLACE_METHOD_HPP
  
