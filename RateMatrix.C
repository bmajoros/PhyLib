/****************************************************************
 RateMatrix.H
 william.majoros@duke.edu

 This is open-source software,
 governed by the Gnu General Public License (GPL) version 3 (see www.opensource.org).
 ****************************************************************/
#include "RateMatrix.H"
#include "BOOM/GSL/LabeledMatrixLoader.H"
#include "BOOM/GSL/Vector.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/GapPatternAlphabet.H"
#include "BOOM/AminoAlphabet.H"
#include "BOOM/Random.H"
#include "JukesCantor.H"
#include "FEL.H"
#include "Kimura2Param.H"
#include "HKY.H"
#include "REV.H"
#include "GapRateMatrix.H"
#include "HalpernBruno.H"
#include <fstream>
using namespace BOOM;


BOOM::Array1D<double> RateMatrix::omitted;


RateMatrix::RateMatrix(MolecularSequenceType seqType,
		       RateMatrixType matrixType)
  : alphabet(seqType==DNA ? (Alphabet&) PureDnaAlphabet::global() : 
	                    (Alphabet&) AminoAlphabet::global()),
    seqType(seqType),
    matrixType(matrixType)
{
  // Allocate matrices and arrays based on alphabet size
  int n=alphabet.getNumElements();
  M.resize(n,n);

  if(seqType==DNA)
    alphabetMap=DropGapMapping(DnaDashDotAlphabet::global(),
			       PureDnaAlphabet::global());
  else
    alphabetMap=AlphabetIdentityMap(alphabet);
}



RateMatrix::RateMatrix(MolecularSequenceType seqType,Alphabet &alphabet,
		       RateMatrixType matrixType)
  : alphabet(alphabet),
    seqType(seqType),
    matrixType(matrixType)
{
  // Allocate matrices and arrays based on alphabet size
  int n=alphabet.getNumElements();
  M.resize(n,n);

  if(seqType==DNA)
    alphabetMap=DropGapMapping(DnaDashDotAlphabet::global(),
			       PureDnaAlphabet::global());
  else
    alphabetMap=AlphabetIdentityMap(alphabet);
}



double RateMatrix::operator()(Symbol a,Symbol b) const
{
  return M(a,b);//alphabetMap(a),alphabetMap(b));
}



double &RateMatrix::operator()(Symbol a,Symbol b)
{
  return M(a,b);//alphabetMap(a),alphabetMap(b));
}



void RateMatrix::printOn(ostream &os) const
{
  int N=M.getNumRows();

  // Determine the spacing
  int maxLen=0;
  for(Symbol i=0 ; i<N ; ++i)
    for(Symbol j=0 ; j<N ; ++j)
      {
	String s=M(i,j);
	int len=s.length();
	if(len>maxLen) maxLen=len;
      }
  int colWidth=maxLen+1;

  // Print the header
  for(int i=0 ; i<colWidth ; ++i) os << ' ';
  for(Symbol i=0 ; i<N ; ++i)
    {
      os << alphabet.lookup(i); //(alphabetMap.unmap(i));
      for(int i=1 ; i<colWidth ; ++i) os << ' ';
    }
  os << endl;

  // Print each row of the matrix
  for(Symbol i=0 ; i<N ; ++i)
    {
      os << alphabet.lookup(i)<<' ';;//alphabetMap.unmap(i)) << ' ';
      for(Symbol j=0 ; j<N ; ++j)
	{
	  os << M(i,j);
	  String s=M(i,j);
	  int len=s.length();
	  int extra=colWidth-len;
	  for(int i=0 ; i<extra ; ++i) os << ' ';
	}
      os << endl;
    }
}



ostream &operator<<(ostream &os,const RateMatrix &M)
{
  M.printOn(os);
  return os;
}



SubstitutionMatrix *RateMatrix::instantiate(double branchLength)
{
  SubstitutionMatrix *s=new SubstitutionMatrix(alphabet,alphabetMap);
  instantiate(branchLength,*s);
  return s;
}



void RateMatrix::instantiate(double branchLength,SubstitutionMatrix &s)
{
  rateMatrixToSubstMatrix(M,branchLength,s.peek());
}



void rateMatrixToSubstMatrix(const GSL::Matrix &rateMatrix,
			     double branchLength,
			     GSL::Matrix &substMatrix)
{
  // Decompose into eigenvectors and eigenvalues
  GSL::Matrix G; // eigenvectors
  GSL::Vector eigenvalues;
  GSL::Vector imagParts; // should be all zeros for practical rate matrices
  rateMatrix.getEigenVectors(G,eigenvalues,imagParts);
  int n=rateMatrix.getNumRows();

  // Prepare Ginverse
  GSL::Matrix Ginverse;
  G.invert(Ginverse);

  // Prepare matrix U (which consists of exp(t*lambda) diagonal values)
  GSL::Matrix U(n,n);
  U.setAllTo(0.0);
  for(int i=0 ; i<n ; ++i)
    U(i,i)=exp(branchLength*eigenvalues[i]);

  // Form the matrix product P(t) = G * U * Ginverse
  GSL::Matrix GU;
  G.times(U,GU);
  GU.times(Ginverse,substMatrix);
}



void RateMatrix::installDiagonals()
{
  int n=M.getNumRows();
  for(int i=0 ; i<n ; ++i)
    {
      double sum=0;
      for(int j=0 ; j<n ; ++j)
	if(j!=i) sum+=M(i,j);
      M(i,i)=-sum;
    }
}



int RateMatrix::getNumRows() const
{
  return M.getNumRows();
}



const Alphabet &RateMatrix::getAlphabet() const
{
  return alphabet;
}



int RateMatrix::numParameters() const
{
  return 0;
}



void RateMatrix::partialDerivative(int parameter,RateMatrix &) const
{
  throw "can't differentiate a non-parametric rate matrix!";
}



void RateMatrix::initMatrix_ACGT(double *X)
{
  RateMatrix &M=*this;
  Symbol A=alphabet.lookup('A');
  Symbol C=alphabet.lookup('C');
  Symbol G=alphabet.lookup('G');
  Symbol T=alphabet.lookup('T');

  M(A,A)=X[0];   M(A,C)=X[1];   M(A,G)=X[2];   M(A,T)=X[3];
  M(C,A)=X[4];   M(C,C)=X[5];   M(C,G)=X[6];   M(C,T)=X[7];
  M(G,A)=X[8];   M(G,C)=X[9];   M(G,G)=X[10];  M(G,T)=X[11];
  M(T,A)=X[12];  M(T,C)=X[13];  M(T,G)=X[14];  M(T,T)=X[15];
}



void RateMatrix::save(const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  save(os);
}



RateMatrix *RateMatrix::load(const BOOM::String &filename)
{
  ifstream is(filename.c_str());
  return load(is);
}



RateMatrix *RateMatrix::load(istream &is)
{
  String line;
  is>>line;
  RateMatrixType T=rateMatrixTypeFromString(line);
  switch(T)
    {
    case MT_JC:   return new JukesCantor(is);
    case MT_FEL:  return new FEL(is);
    case MT_KIM:  return new Kimura2Param(is);
    case MT_HKY:  return new HKY(is);
    case MT_REV:  return new REV(is);
    case MT_GAP:  return new GapRateMatrix(is);
    case MT_HB:   return new HalpernBruno(is);
    }
  throw String("unknown rate matrix type")+T;
}



RateMatrix *RateMatrix::random(MolecularSequenceType seqType,
			       RateMatrixType matrixType,
			       BOOM::Array1D<double> equilibriumFreqs)
{
  if(seqType!=DNA) throw "RateMatrix::Random(PROTEIN) not implemented yet";

  if(matrixType==MT_GAP) return GapRateMatrix::random();

  // Generate random equilib freqs if the caller didn't provide them
  if(equilibriumFreqs.size()==0)
    {
      int n;
      if(seqType==DNA)
	if(matrixType==MT_GAP)
	  n=GapPatternAlphabet::global().size();
	else
	  n=PureDnaAlphabet::global().size();
      else
	n=AminoAlphabet::global().size();
      equilibriumFreqs.resize(n);
      for(int i=0 ; i<n ; ++i)
	equilibriumFreqs[i]=RandomFloat(0.2,0.3); // 20% - 30%
      double total=0.0;
      for(int i=0 ; i<n ; ++i)
	total+=equilibriumFreqs[i];
      for(int i=0 ; i<n ; ++i)
	equilibriumFreqs[i]/=total;
    }

  // Generate model parameters
  double transitions=  RandomFloat(0.1,0.5);
  double transversions=RandomFloat(0.1,0.3);
  double alpha=        RandomFloat(0.1,0.8);
  double beta=         RandomFloat(0.1,0.3);
  double kappa=        RandomFloat(0.1,0.3);
  double chi=          RandomFloat(0.1,0.3);
  double omega=        RandomFloat(0.1,0.8);
  double tau=          RandomFloat(0.1,0.3);
  double mutationRate= RandomFloat(0.1,0.3);

  // Build the matrix
  switch(matrixType)
    {
    case MT_JC: return new JukesCantor(mutationRate);
    case MT_FEL:return new FEL(equilibriumFreqs,mutationRate);
    case MT_KIM:return new Kimura2Param(transitions,transversions);
    case MT_HKY:return new HKY(equilibriumFreqs,transitions,transversions);
    case MT_REV:return new REV(equilibriumFreqs,alpha,beta,kappa,chi,
			       omega,tau);
    default: throw "RateMatrix::random";
    }
}



RateMatrixType RateMatrix::getMatrixType()
{
  return matrixType;
}



void RateMatrix::rescale(double factor)
{
  M.scale(factor);
  /*
  int n=M.getNumRows();
  for(int i=0 ; i<n ; ++i)
    for(int j=0 ; j<n ; ++j)
      M(i,j)*=factor;
  */
}



const GSL::Matrix &RateMatrix::peek() const
{
  return M;
}



RateMatrix *RateMatrix::complement() const
{
  const RateMatrix &self=*this;
  RateMatrix &other=*clone();//new RateMatrix(seqType,alphabet,matrixType);
  int n=alphabet.getNumElements();
  for(Symbol a=0 ; a<n ; ++a) {
    Symbol aC=alphabet.complement(a);
    for(Symbol b=0 ; b<n ; ++b) {
      Symbol bC=alphabet.complement(b);
      other(a,b)=self(aC,bC);
    }
  }
  return &other;
}


