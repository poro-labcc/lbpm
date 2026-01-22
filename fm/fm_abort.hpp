//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
//
// Function to abort the program when an error occours.
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef ABORT_HPP
#define ABORT_HPP


#include <iostream>
#include <cstdlib>
using namespace std;


#define abort_fm( error ){                                         \
  cerr << "!!  File: " << __FILE__ << "\n"                   \
       << "!!  Line  : " << __LINE__ << "\n"                   \
       << "!!  " << error << "\n\a!!  ABORTING ..." << endl;    \
  exit(1);                                                      \
  }


#endif /* ABORT_HPP */
