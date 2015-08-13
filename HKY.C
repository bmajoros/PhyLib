/****************************************************************
 HKY.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "HKY.H"
#include "BOOM/DnaAlphabet.H"

HKY::HKY(BOOM::Array1D<double> eq,double alpha,double beta)
  : RateMatrix(DNA,MT_HKY),
    eq(eq),
    alpha(alpha),
    beta(beta)
{
  // ctor

  init();
}



RateMatrix *HKY::clone() const {
    return new HKY(eq,alpha,beta);
}



HKY::HKY(istream &is)
  : RateMatrix(DNA)
{
  is>>alpha>>beta>>eq;
  init();
}



void HKY::save(ostream &os)
{
  os<<"HKY"<<endl;
  os<<alpha<<" "<<beta<<" "<<eq<<endl;
}



void HKY::init()
{
  // Requirement: alpha+2*beta<1.  This block of code scales the
  // values so that this constrain holds, without changing the
  // ratio transitionP/transversionP:
  double alphaPlusTwoBeta=(alpha+2*beta)/0.999;
  alpha/=alphaPlusTwoBeta;
  beta/=alphaPlusTwoBeta;

  RateMatrix &self=*this;
  const Alphabet &alphabet=getAlphabet();
  Symbol A=alphabet.lookup('A');
  Symbol C=alphabet.lookup('C');
  Symbol G=alphabet.lookup('G');
  Symbol T=alphabet.lookup('T');
  double piA=eq[A], piC=eq[C], piG=eq[G], piT=eq[T];

  self(A,C)=beta*piC;
  self(A,G)=alpha*piG;
  self(A,T)=beta*piT;

  self(C,A)=beta*piA;
  self(C,G)=beta*piG;
  self(C,T)=alpha*piT;

  self(G,A)=alpha*piA;
  self(G,C)=beta*piC;
  self(G,T)=beta*piT;

  self(T,A)=beta*piA;
  self(T,C)=alpha*piC;
  self(T,G)=beta*piG;

  installDiagonals();
}



int HKY::numParameters() const
{
  return 2;
}



void HKY::partialDerivative(int parameter,RateMatrix &D) const
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

  switch(parameter)
    {
    case 0: // alpha
      {
	double array[16]=
	  {
	    -piG,     0,   piG,     0,
	       0,  -piT,     0,   piT,
	     piA,     0,  -piA,     0,
	       0,   piC,     0,  -piC
	  };
	initMatrix_ACGT(array);
      }
      break;
    case 1: // beta
      {
	double array[16]=
	  {
	    -piC-piT,      piC,        0,       piT,
	         piA, -piA-piG,      piG,         0,
	           0,      piC, -piC-piT,       piT,
	         piA,        0,      piG,  -piA-piG
	  };
	initMatrix_ACGT(array);
      }
      break;
    }
}



RateMatrixType HKY::getType()
{
  return MT_HKY;
}



double HKY::getIthParm(int i) const
{
  return i==0 ? alpha : beta;
}



void HKY::setIthParm(int i,double parm)
{
  if(i==0) alpha=parm; else beta=parm;
}



RateMatrix *HKY::average(const Array1D<RateMatrix*> &matrices) const {
    int n=matrices.size();
    if(n==0) return NULL;
    int eqDim=dynamic_cast<HKY*>(matrices[0])->eq.size();
    BOOM::Array1D<double> eq(eqDim);
    eq.setAllTo(0.0);
    double alpha=0.0, beta=0.0;
    for(int i=0 ; i<n ; ++i) {
        HKY *M=dynamic_cast<HKY*>(matrices[i]);
        alpha+=M->alpha;
        beta+=M->beta;
        for(int j=0 ; j<eqDim ; ++j) eq[j]+=M->eq[j];
    }
    alpha/=n;
    beta/=n;
    for(int j=0 ; j<eqDim ; ++j) eq[j]/=n;
    return new HKY(eq,alpha,beta);
}



void HKY::addNoise(GSL::ContinuousDistribution &d) {
    alpha+=d.random();
    beta+=d.random();

    if(alpha<=0) alpha=0.01;
    else if(alpha>=1) alpha=0.99;

    if(beta<=0) beta=0.01;
    else if(beta>=1) beta=0.99;

    init();
}



