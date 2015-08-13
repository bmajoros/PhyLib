/****************************************************************
 FEL.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "FEL.H"
#include "BOOM/DnaAlphabet.H"
using namespace std;
using namespace BOOM;

FEL::FEL(BOOM::Array1D<double> eq,double alpha)
  : RateMatrix(DNA,MT_FEL),
    eq(eq),
    alpha(alpha)
{
  // ctor

  init();
}



RateMatrix *FEL::clone() const {
    return new FEL(eq,alpha);
}



FEL::FEL(istream &is)
  : RateMatrix(DNA)
{
  is>>alpha>>eq;
  init();
}



void FEL::save(ostream &os)
{
  os<<"FEL"<<endl;
  os<<alpha<<" "<<eq<<endl;
}



void FEL::init()
{
  RateMatrix &self=*this;
  const Alphabet &alphabet=getAlphabet();
  int n=alphabet.size();
  for(Symbol to=0 ; to<n ; ++to)
    {
      double entry=alpha*eq[to];
      for(Symbol from=0 ; from<n ; ++from)
	if(from!=to)
	  self(from,to)=entry;
    }
  installDiagonals();
}



int FEL::numParameters() const
{
  return 1;
}



void FEL::partialDerivative(int parameter,RateMatrix &D) const
{
  Alphabet &alpha=DnaAlphabet::global();
  Symbol A=alpha.lookup('A');
  Symbol C=alpha.lookup('C');
  Symbol G=alpha.lookup('G');
  Symbol T=alpha.lookup('T');
  double piA=eq[A];
  double piC=eq[C];
  double piG=eq[G];
  double piT=eq[T];
  double array[16]=
    {
      -piC-piG-piT,      piC,          piG,          piT,
           piA,     -piA-piG-piT,      piG,          piT,
           piA,          piC,     -piA-piC-piT,      piT,
           piA,          piC,          piG,     -piA-piC-piG
    };
  initMatrix_ACGT(array);
}



RateMatrixType FEL::getType()
{
  return MT_FEL;
}



double FEL::getIthParm(int i) const
{
  return alpha;
}



void FEL::setIthParm(int i,double parm)
{
  alpha=parm;
  init();
}



RateMatrix *FEL::average(const Array1D<RateMatrix*> &matrices) const {
    int n=matrices.size();
    if(n==0) return NULL;
    int eqDim=dynamic_cast<FEL*>(matrices[0])->eq.size();
    BOOM::Array1D<double> eq(eqDim);
    eq.setAllTo(0.0);
    double alpha=0.0;
    for(int i=0 ; i<n ; ++i) {
        FEL *M=dynamic_cast<FEL*>(matrices[i]);
        alpha+=M->alpha;
        for(int j=0 ; j<eqDim ; ++j) eq[j]+=M->eq[j];
    }
    alpha/=n;
    for(int j=0 ; j<eqDim ; ++j) eq[j]/=n;
    return new FEL(eq,alpha);
}



void FEL::addNoise(GSL::ContinuousDistribution &d) {
    alpha+=d.random();
    if(alpha<=0) alpha=0.01;
    else if(alpha>=1) alpha=0.99;

    init();
}



