/****************************************************************
 JukesCantor.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "JukesCantor.H"
using namespace BOOM;


JukesCantor::JukesCantor(double alpha)
  : RateMatrix(DNA,MT_JC), alpha(alpha)
{
  // ctor

  init();
}



RateMatrix *JukesCantor::clone() const {
    return new JukesCantor(alpha);
}



void JukesCantor::init()
{
  RateMatrix &self=*this;
  const Alphabet &alphabet=getAlphabet();
  int n=alphabet.size();
  for(Symbol i=0 ; i<n ; ++i)
    for(Symbol j=0 ; j<n ; ++j)
      if(j!=i) 
	self(i,j)=alpha;
  
  installDiagonals();
}



JukesCantor::JukesCantor(istream &is)
  : RateMatrix(DNA)
{
  is>>alpha;
  init();
}



void JukesCantor::save(ostream &os)
{
  os<<"JC"<<endl;
  os<<alpha<<endl;
}



int JukesCantor::numParameters() const
{
  return 1;
}



void JukesCantor::partialDerivative(int parameter,RateMatrix &D) const
{
  double M[16]=
  {
    -3,   1,   1,   1,
     1,  -3,   1,   1,
     1,   1,  -3,   1,
     1,   1,   1,  -3
  };
  D.initMatrix_ACGT(M);
}



RateMatrixType JukesCantor::getType()
{
  return MT_JC;
}



double JukesCantor::getIthParm(int i) const
{
  return alpha;
}



void JukesCantor::setIthParm(int i,double parm)
{
  alpha=parm;
  init();
}



RateMatrix *JukesCantor::average(const Array1D<RateMatrix*> &matrices) const {
    int n=matrices.size();
    double alpha=0;
    for(int i=0 ; i<n ; ++i)
        alpha+=dynamic_cast<JukesCantor*>(matrices[i])->alpha;
    alpha/=n;
    return new JukesCantor(alpha);
}



void JukesCantor::addNoise(GSL::ContinuousDistribution &d) {
    alpha+=d.random();
    if(alpha<=0) alpha=0.01;
    else if(alpha>=1) alpha=0.99;

    init();
}


