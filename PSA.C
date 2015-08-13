/****************************************************************
 PSA.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "PSA.H"
using namespace std;
using namespace BOOM;


PSA::PSA(const String &filename)
  : file(filename,"r")
{
}



int PSA::getLength()
{
  return file.getSize()/sizeof(double);
}



double PSA::getLogLikelihood(long begin,long end) // half-open interval: [b,e)
{
  double beginValue;
  if(begin>0) {
    file.seek(begin-1);
    beginValue=file.readDouble();
  }
  else beginValue=0.0;
  file.seek(end-1);
  double endValue=file.readDouble();
  return endValue-beginValue;
}



void PSA::close()
{
  file.close();
}


