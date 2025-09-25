//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
//
// Function to abort the program when an error occours.
//________________________________________________________
//A.Z. - 03/05 => Creation
//       11/05 => Uses #define instead of a function.
//                Shows the file and line where the error occoured.
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef ABORTA_HPP
#define ABORTA_HPP


#include <iostream>
#include <cstdlib>
using namespace std;


#define aborta( erro ){                                         \
  cerr << "!!  File: " << __FILE__ << "\n"                   \
       << "!!  Line  : " << __LINE__ << "\n"                   \
       << "!!  " << erro << "\n\a!!  ABORTING ..." << endl;    \
  exit(1);                                                      \
  }


#endif /* ABORTA_HPP */
