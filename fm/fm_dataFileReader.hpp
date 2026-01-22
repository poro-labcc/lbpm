//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// Class DataFileReader
//
// This class provides a convenient interface for reading
// structured data files. It automatically skips comments
// and keeps track of the current line number, column number,
// and the total number of columns in each line.

//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|

// THIS FILE IS NOT NECESSARY ANYMORE.



#ifndef DATAFILEREADER_H
#define DATAFILEREADER_H

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "fm_abort.hpp"
// #include "fm_numbers.hpp"
#include "fm_types.hpp"
using namespace std;


class DataFileReader{
  template<class T>
  friend DataFileReader &operator >> ( DataFileReader &, T & );

public:
  DataFileReader( MTcs &, MTcs & = "#" );
  ~DataFileReader(){ this->close(); }

public:
  void   commentarySymbol( MTcs &cs )       { _commentarySymbol = cs;   }
  string commentarySymbol(   void   ) const { return _commentarySymbol; }
  int    lineNumber      (   void   ) const { return _lineNumber;       }
  int    columnNumber    (   void   ) const { return _columnNumber;     }
  int    numberOfColumns (   void   ) const { return _numberOfColumns;  }

  void open( MTcs &,  MTcs & = "#");

  template<class T> T readColumn ( MTci & );

  void   close(void){  if( _FILE.is_open() ) _FILE.close(); }

  bool operator++();     

  const DataFileReader &operator += ( MTci &n ){
    for( int i=0; i<n; i++ ) ++(*this);
    return (*this); 
  }

  void rewind(void);

private:
  int      _lineNumber;
  int      _columnNumber;
  int      _numberOfColumns;
  string   _commentarySymbol;
  string   _filename;
  MTvs     _lineContent;
  ifstream _FILE;
};





//------------------------------------------------------------------------------
//   Constructor
//
DataFileReader::DataFileReader( MTcs &file, MTcs &cs ){
  this->open(file,cs);
}





//------------------------------------------------------------------------------
// Opens a file and initializes the reader state.
//
void DataFileReader::open( MTcs &file, MTcs &cs ){
  _commentarySymbol = cs;
  _filename=file;
  this->rewind();
}



//------------------------------------------------------------------------------
// Reads the value from column C (1-based index).
// Aborts if the column does not exist.
//
template<class T>
T DataFileReader::readColumn( MTci &C ){
  if(C<1){
    abort_fm("Columns are indexed starting from 1.");
  }
  if(C>_numberOfColumns){
    string msg="Line "+ntos(_lineNumber)+" has no "+ntos(C)+" column.";
    msg+=" Valid columns range from 1 to ";
    msg+=ntos(_numberOfColumns) + ".";
    abort_fm(msg);
  }
  T v;
  stringstream ss( _lineContent[C-1] );
  ss >> v;
  return v;
}


//------------------------------------------------------------------------------
// Overloaded extraction operator.
// 
template<class T>
DataFileReader &operator >> ( DataFileReader &s, T &v ){

  v = s.template readColumn<T>(s._columnNumber);

  s._columnNumber++;
  if( s._columnNumber == (s._numberOfColumns+1) ) s._columnNumber = 1;

  return s;
}



//------------------------------------------------------------------------------
// Pre-increment operator.
//
bool DataFileReader::operator++(){

  string saux1,saux2;
  _columnNumber = 1;
  _lineContent  = MTvs();
  while(  _lineContent.empty() && !_FILE.eof() ){

    _lineNumber++;
    string completeLine;
    getline( _FILE, completeLine );


    unsigned int p = completeLine.find(_commentarySymbol);


    if( p != string::npos )
      saux1.assign( completeLine, 0, p );
    else                 
      saux1.assign( completeLine );


    stringstream ssaux(saux1);
    while( ssaux >> saux2 ) _lineContent.push_back(saux2);
    _numberOfColumns = _lineContent.size();
  }

  return (_FILE.eof() && _lineContent.size()==0);
}



//------------------------------------------------------------------------------
//   Rewinds the reader to the beginning of the file.
//
void DataFileReader::rewind(void){
  if( _FILE.is_open() ) this->close();  
  _FILE.open( _filename.c_str() );      
  if( !_FILE.is_open() )  abort_fm("Unable to open file "+_filename);

  _lineNumber      = 0; 
  _columnNumber    = 1;
  _numberOfColumns = 0;
 }



#endif // DATAFILEREADER_H
