/****************************************************************
 HalpernBruno.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "HalpernBruno.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;



HalpernBruno::HalpernBruno(RateMatrix &backgroundMutationRate,
			   const Array1D<double> &eqFreqs)
  : eq(eqFreqs), Q(backgroundMutationRate.clone()), RateMatrix(DNA,MT_HB)
{
  // ctor

  init();
}



HalpernBruno::HalpernBruno(istream &is)
  : RateMatrix(DNA,MT_HB)
{
  is>>eq;
  Q=RateMatrix::load(is);
  init();
}



HalpernBruno::~HalpernBruno()
{
  delete Q;
}



void HalpernBruno::save(ostream &os)
{
  os<<"HB"<<endl;
  os<<eq<<endl;
  Q->save(os);
}



double HalpernBruno::getIthParm(int i) const
{
  return Q->getIthParm(i);
}



int HalpernBruno::numParameters() const
{
  return Q->numParameters();
}



void HalpernBruno::setIthParm(int i,double parm)
{
  Q->setIthParm(i,parm);
}



RateMatrix *HalpernBruno::average(const Array1D<RateMatrix*> &) const
{
  throw "Halpern Bruno::average() -- not implemented";
}



void HalpernBruno::addNoise(GSL::ContinuousDistribution &)
{
  throw "HalpernBruno::addNoise() -- not implemented";
}



RateMatrix *HalpernBruno::clone() const
{
  return new HalpernBruno(*Q,eq);
}



void HalpernBruno::init()
{
  for(int a=0 ; a<4 ; ++a) {
    double sum=0;
    for(int b=0 ; b<4 ; ++b) {
      if(b==a) continue;
      double Qab=(*Q)(a,b), Qba=(*Q)(b,a);
      double ratio=eq[b]*Qba/(eq[a]*Qab);
      M(a,b)=(ratio==1.0 ? Qab : Qab*log(ratio)/(1-1/ratio));
      //if(!isFinite(M(a,b))) INTERNAL_ERROR;
      sum+=M(a,b);
    }
    M(a,a)=-sum;
    //if(!isFinite(sum)) INTERNAL_ERROR;
  }
}


