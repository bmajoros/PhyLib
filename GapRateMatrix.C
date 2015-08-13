/****************************************************************
 GapRateMatrix.H
 william.majoros@duke.edu

 This is open-source software,
 governed by the Gnu General Public License (GPL) version 3 (see www.opensource.org).
 ****************************************************************/
#include "GapRateMatrix.H"
#include "BOOM/GSL/LabeledMatrixLoader.H"
#include "BOOM/GSL/Vector.H"
#include "BOOM/Random.H"
#include <fstream>
using namespace BOOM;



GapRateMatrix::GapRateMatrix()
  : RateMatrix(DNA,GapPatternAlphabet::global(),MT_GAP)
{
  alphabetMap=GapPatternAlphaMap::global();
}



GapRateMatrix::GapRateMatrix(istream &is)
  : RateMatrix(DNA,GapPatternAlphabet::global(),MT_GAP)
{
  alphabetMap=GapPatternAlphaMap::global();
  load(is);
  init();
}



int GapRateMatrix::numParameters() const
{
  return 2;
}



void GapRateMatrix::partialDerivative(int parameter,GapRateMatrix &) const
{
  throw "not implemented yet -- sorry about that";
}



GapRateMatrix *GapRateMatrix::random()
{
  GapRateMatrix *m=new GapRateMatrix;
  m->dashParm=RandomFloat(0.05,0.5);
  m->dotParm= RandomFloat(0.05,0.5);
  m->init();
  return m;
}



void GapRateMatrix::save(ostream &os)
{
  os<<"GAP"<<endl;
  os<<dashParm<<"\n"<<dotParm<<endl;
}



void GapRateMatrix::load(istream &is)
{
  String type;
  is>>type;
  if(type!="GAP") 
    throw String("Trying to load a rate matrix of type ")+type+
      " into a GapRateMatrix";
  is>>dashParm>>dotParm;
}



double GapRateMatrix::getIthParm(int i) const
{
  switch(i)
    {
    case 0: return dashParm;
    case 1: return dotParm;
    }
}



void GapRateMatrix::setIthParm(int i,double parm)
{
  switch(i)
    {
    case 0: dashParm=parm; break;
    case 1: dotParm=parm; break;
    }
}



RateMatrix *GapRateMatrix::average(const Array1D<RateMatrix*> &a) const
{
  throw "GapRateMatrix::average() not implemented";
}



void GapRateMatrix::addNoise(GSL::ContinuousDistribution &d)
{
  throw "GapRateMatrix::addNoise() not implemented";
}



GapRateMatrix *GapRateMatrix::clone() const
{
  GapRateMatrix *other=new GapRateMatrix;
  other->dashParm=dashParm;
  other->dotParm=dotParm;
  return other;
}



void GapRateMatrix::init()
{
  RateMatrix &self=*this;
  Symbol N=alphabet.lookup('N');
  Symbol dash=alphabet.lookup('-');
  Symbol dot=alphabet.lookup('.');
  M(N,dash)=dashParm;
  M(N,dot)=dotParm;
  M(dash,N)=0;
  M(dash,dot)=0;
  M(dot,N)=0;
  M(dot,dash)=0;

  installDiagonals();  
}



void GapRateMatrix::printOn(ostream &os) const
{
  os<<"P(-)="<<dashParm<<", P(.)="<<dotParm<<endl;
}



ostream &operator<<(ostream &os,const GapRateMatrix &M)
{
  M.printOn(os);
  return os;
}




