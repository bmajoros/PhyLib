/****************************************************************
 REV.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "REV.H"
#include "BOOM/DnaAlphabet.H"


REV::REV(BOOM::Array1D<double> eq,double alpha,double beta,double kappa,
	 double chi,double omega,double tau)
  : RateMatrix(DNA,MT_REV),
    eq(eq),
    alpha(alpha), beta(beta), kappa(kappa), chi(chi), omega(omega), tau(tau)
{
  // ctor

  init();
}



RateMatrix *REV::clone() const {
    return new REV(eq,alpha,beta,kappa,chi,omega,tau);
}



REV::REV(istream &is)
  : RateMatrix(DNA)
{
  is>>alpha>>beta>>kappa>>chi>>omega>>tau>>eq;
  init();
}



void REV::save(ostream &os)
{
  os<<"REV"<<endl;
  os<<alpha<<" "<<beta<<" "<<kappa<<" "<<chi<<" "<<omega<<" "
    <<tau<<"\n"<<eq<<endl;
}



void REV::init()
{
  RateMatrix &self=*this;
  const Alphabet &alphabet=getAlphabet();
  Symbol A=alphabet.lookup('A');
  Symbol C=alphabet.lookup('C');
  Symbol G=alphabet.lookup('G');
  Symbol T=alphabet.lookup('T');
  double piA=eq[A], piC=eq[C], piG=eq[G], piT=eq[T];

  self(A,C)=beta*piC;
  self(A,G)=alpha*piG;
  self(A,T)=chi*piT;

  self(C,A)=beta*piA;
  self(C,G)=kappa*piG;
  self(C,T)=omega*piT;

  self(G,A)=alpha*piA;
  self(G,C)=kappa*piC;
  self(G,T)=tau*piT;

  self(T,A)=chi*piA;
  self(T,C)=omega*piC;
  self(T,G)=tau*piG;

  installDiagonals();
}



int REV::numParameters() const
{
  return 6;
}



void REV::partialDerivative(int parameter,RateMatrix &D) const
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
    case 0:
      {
	double array[16]=
	  {
	    -piG,          0,        piG,          0,
	       0,          0,          0,          0,
	     piA,          0,       -piA,          0,
	       0,          0,          0,          0
	  };
	initMatrix_ACGT(array);
      }
      break;

    case 1:
      {
	double array[16]=
	  {
	    -piC,        piC,          0,          0,
	     piA,       -piA,          0,          0,
	       0,          0,          0,          0,
	       0,          0,          0,          0
	  };
	initMatrix_ACGT(array);
      }
      break;

    case 2:
      {
	double array[16]=
	  {
	       0,          0,          0,          0,
	       0,       -piG,        piG,          0,
	       0,        piC,       -piC,          0,
	       0,          0,          0,          0
	  };
	initMatrix_ACGT(array);
      }
      break;

    case 3:
      {
	double array[16]=
	  {
	    -piT,          0,          0,        piT,
	       0,          0,          0,          0,
	       0,          0,          0,          0,
	     piA,          0,          0,       -piA
	  };
	initMatrix_ACGT(array);
      }
      break;

    case 4:
      {
	double array[16]=
	  {
	       0,          0,          0,          0,
	       0,       -piT,          0,        piT,
	       0,          0,          0,          0,
	       0,        piC,          0,       -piC
	  };
	initMatrix_ACGT(array);
      }
      break;

    case 5:
      {
	double array[16]=
	  {
	       0,          0,          0,          0,
	       0,          0,          0,          0,
	       0,          0,       -piT,        piT,
	       0,          0,        piG,       -piG
	  };
	initMatrix_ACGT(array);
      }
      break;
    }
}



RateMatrixType REV::getType()
{
  return MT_REV;
}



double REV::getIthParm(int i) const
{
  switch(i)
    {
    case 0: return alpha;
    case 1: return beta;
    case 2: return kappa;
    case 3: return chi;
    case 4: return omega;
    case 5: return tau;
    }
}



void REV::setIthParm(int i,double parm)
{
  switch(i)
    {
    case 0: alpha=parm; break;
    case 1: beta=parm;  break;
    case 2: kappa=parm; break;
    case 3: chi=parm;   break;
    case 4: omega=parm; break;
    case 5: tau=parm;   break;
    }
  //init();
}



//  double alpha, beta, kappa, chi, omega, tau;

RateMatrix *REV::average(const Array1D<RateMatrix*> &matrices) const {
    int n=matrices.size();
    if(n==0) return NULL;
    int eqDim=dynamic_cast<REV*>(matrices[0])->eq.size();
    BOOM::Array1D<double> eq(eqDim);
    eq.setAllTo(0.0);
    double alpha=0,beta=0,kappa=0,chi=0,omega=0,tau=0;
    for(int i=0 ; i<n ; ++i) {
        REV *M=dynamic_cast<REV*>(matrices[i]);
        alpha+=M->alpha;
        beta+=M->beta;
        kappa+=M->kappa;
        chi+=M->chi;
        omega+=M->omega;
        tau+=M->tau;
        for(int j=0 ; j<eqDim ; ++j) eq[j]+=M->eq[j];
    }
    alpha/=n;
    beta/=n;
    kappa/=n;
    chi/=n;
    omega/=n;
    tau/=n;
    for(int j=0 ; j<eqDim ; ++j) eq[j]/=n;
    return new REV(eq,alpha,beta,kappa,chi,omega,tau);
}



void REV::bound(double &param) {
    if(param<=0) param=0.01;
    else if(param>=1) param=0.99;
}



void REV::sumToLessThanOne(double &a,double &b,double &c) {
    double sum=a+b+c;
    if(sum>=1) {
        double d=(sum-1)*0.36; // slightly more than 1/3
        a-=d; if(a<=0) a=0.01;
        b-=d; if(b<=0) b=0.01;
        c-=d; if(c<=0) c=0.01;
    }
}



void REV::addNoise(GSL::ContinuousDistribution &d) {
    alpha+=d.random(); bound(alpha);
    beta+=d.random();  bound(beta);
    kappa+=d.random(); bound(kappa);
    chi+=d.random();   bound(chi);
    omega+=d.random(); bound(omega);
    tau+=d.random();   bound(tau);

    sumToLessThanOne(alpha,beta,chi);
    sumToLessThanOne(beta,kappa,omega);
    sumToLessThanOne(alpha,kappa,tau);
    sumToLessThanOne(chi,omega,tau);
    
    init();
}

