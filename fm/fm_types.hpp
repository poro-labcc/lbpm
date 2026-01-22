//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
//
//   Most used types
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef FMTYPES_H
#define FMTYPES_H

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

#endif /* FMTYPES_H */
