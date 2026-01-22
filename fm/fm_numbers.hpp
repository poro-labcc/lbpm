//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
//
// Library with many functions to manipulate numbers.
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|

// THIS FILE IS NOT NECESSARY ANYMORE.


#ifndef NUMEROS_H
#define NUMEROS_H


#include <deque>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <iostream>
using namespace std;


template < class N > inline string ntos(const N &);
template < class N > inline string ntos(const N &, const N &);
void binario(int,int,deque<int> &);




//==============================================================================
// Converts a number into a string
template < class N >
inline string ntos(const N &n){
  ostringstream s;
  s << n;
  return s.str();
}




//==============================================================================
// Given a number n, formats it with the number of figures of N:
// Example ( N=5000 -> 4 figures )
// 1    --> 0001
// 25   --> 0025
// 125  --> 0125
// 3125 --> 3125
template < class N >
inline string ntos(const N&model, const N &n){

  // Number of model number figures
  stringstream smodel;
  smodel << model;
  const int length = smodel.str().length();

  // Formats output
  stringstream sn;
  sn.width(length);
  sn.fill('0');

  // Generates number
  sn << n;
  return sn.str();
}




//==============================================================================
// Given a number n, formats it to binary with N figures
// Example: N = 5
// n=0 --> 0 0 0 0 0
// n=1 --> 0 0 0 0 1
// n=2 --> 0 0 0 1 0
// void binario(int n, int N, deque<int> &bin){
//   bin.clear(); 

//   while( n>0 ){        
//     int mod = n%2;
//     n -= mod;
//     n /= 2;
//     bin.push_front(mod);
//   }

//   int size=bin.size(); 
//   for( int i=size; i<N; i++ ){
//     bin.push_front(0);
//   }

//   if(size>N){        
//     cerr << "\n\aIt is impossible to write the number in binary with only"
// 	  << N << " figures! ABORTING ..." << endl;
//     exit(1);
//   }
// }



#endif /* NUMEROS_H */
