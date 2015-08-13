/****************************************************************
 PSA.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_PSA_H
#define INCL_PSA_H
#include <iostream>
#include "BOOM/File.H"
using namespace std;
using namespace BOOM;

// Implements a prefix-sum array on disk

class PSA
{
  File file;
public:
  PSA(const String &filename);
  int getLength();
  double getLogLikelihood(long begin,long end); // half-open interval: [b,e)
  void close();
};

#endif

