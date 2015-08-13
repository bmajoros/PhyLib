/****************************************************************
 Kimura2Param.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "Kimura2Param.H"


Kimura2Param::Kimura2Param(double transitionP,double transversionP)
  : RateMatrix(DNA,MT_KIM),
    transitionP(transitionP),
    transversionP(transversionP)
{
  // ctor

  init();
}



RateMatrix *Kimura2Param::clone() const {
    return new Kimura2Param(transitionP,transversionP);
}



Kimura2Param::Kimura2Param(istream &is)
  : RateMatrix(DNA)
{
  is>>transitionP>>transversionP;
  init();
}



void Kimura2Param::save(ostream &os)
{
  os<<"KIM"<<endl;
  os<<transitionP<<" "<<transversionP<<endl;
}



void Kimura2Param::init()
{
  // Requirement: alpha+2*beta<1.  This block of code scales the
  // values so that this constrain holds, without changing the
  // ratio transitionP/transversionP:
  double alphaPlusTwoBeta=(transitionP+2*transversionP)/0.999;
  transitionP/=alphaPlusTwoBeta;
  transversionP/=alphaPlusTwoBeta;

  // Install values into matrix
  RateMatrix &self=*this;
  const Alphabet &alphabet=getAlphabet();
  Symbol A=alphabet.lookup('A');
  Symbol C=alphabet.lookup('C');
  Symbol G=alphabet.lookup('G');
  Symbol T=alphabet.lookup('T');
  self(A,G)=self(G,A)=self(C,T)=self(T,C)=transitionP;
  self(A,C)=self(A,T)=self(C,A)=self(C,G)=self(G,T)=self(G,C)=self(T,A)=
    self(T,G)=transversionP;
  installDiagonals();
}



int Kimura2Param::numParameters() const
{
  return 2;
}



void Kimura2Param::partialDerivative(int parameter,RateMatrix &D) const
{
  switch(parameter)
    {
    case 0: // alpha
      {
	double M[16]=
	  {
	    -1,   0,   1,   0,
	     0,  -1,   0,   1,
	     1,   0,  -1,   0,
	     0,   1,   0,  -1
	  };
	D.initMatrix_ACGT(M);
      }
      break;

    case 1: // beta
      {
	double M[16]=
	  {
	    -2,   1,   0,   1,
	     1,  -2,   1,   0,
	     0,   1,  -2,   1,
	     1,   0,   1,  -2
	  };
	D.initMatrix_ACGT(M);
      }
      break;
    }
}



RateMatrixType Kimura2Param::getType()
{
  return MT_KIM;
}



double Kimura2Param::getIthParm(int i) const
{
  return i==0 ? transitionP : transversionP;
}



void Kimura2Param::setIthParm(int i,double parm)
{
  if(i==0) transitionP=parm; else transversionP=parm;
}


RateMatrix *Kimura2Param::average(const Array1D<RateMatrix*> &matrices) const {
    int n=matrices.size();
    double transitionP=0, transversionP=0;
    for(int i=0 ; i<n ; ++i) {
        Kimura2Param *M=dynamic_cast<Kimura2Param*>(matrices[i]);
        transitionP+=M->transitionP;
        transversionP+=M->transversionP;
    }
    transitionP/=n;
    transversionP/=n;
    return new Kimura2Param(transitionP,transversionP);
}



void Kimura2Param::addNoise(GSL::ContinuousDistribution &d) {
    transitionP+=d.random();
    transversionP+=d.random();

    if(transitionP<=0) transitionP=0.01;
    else if(transitionP>=1) transitionP=0.99;
    if(transversionP<=0) transversionP=0.01;
    else if(transversionP>=1) transversionP=0.99;
    
    init();
}


