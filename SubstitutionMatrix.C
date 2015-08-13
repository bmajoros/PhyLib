/***********************************************************************
 SubstitutionMatrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ***********************************************************************/
#include "SubstitutionMatrix.H"
#include "BOOM/GSL/LabeledMatrixLoader.H"
#include "BOOM/GSL/Vector.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/AminoAlphabet.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/DnaDashDotAlphabet.H"
using namespace BOOM;


SubstitutionMatrix::SubstitutionMatrix(const SubstitutionMatrix &other)
  : alphabet(other.alphabet), alphabetMap(other.alphabetMap),
    M(other.M), seqType(other.seqType), generators(other.generators),
    logSpace(false)
{
  // copy ctor  
}



SubstitutionMatrix::SubstitutionMatrix(const String &filename,
				       MolecularSequenceType seqType)
  : alphabet(seqType==DNA ? (Alphabet&) PureDnaAlphabet::global() : 
	     (Alphabet&) AminoAlphabet::global()),
    logSpace(false)
{
  if(seqType==DNA)
    alphabetMap=DropGapMapping(DnaDashDotAlphabet::global(),
			       PureDnaAlphabet::global());
  else
    alphabetMap=AlphabetIdentityMap(alphabet);

  // Allocate matrices and arrays based on alphabet size
  int n=alphabet.size();//alphabetMap.getRangeSize();
  M.resize(n,n);

  // Load data
  GSL::LabeledMatrixLoader::load(alphabet,alphabetMap,filename,M);
}



SubstitutionMatrix::SubstitutionMatrix(const Alphabet &a,const AlphabetMap &m)
  : alphabet(a), alphabetMap(m), logSpace(false)
{
  // Allocate matrices and arrays based on alphabet size
  int n=alphabet.size();//alphabetMap.getRangeSize();
  M.resize(n,n);
}



bool SubstitutionMatrix::isInLogSpace() const
{
  return logSpace;
}



GSL::Matrix &SubstitutionMatrix::peek()
{
  return M;
}



double SubstitutionMatrix::operator()(Symbol a,Symbol b) 
  const
{
  return M(a,b);//alphabetMap(a),alphabetMap(b));
}



double &SubstitutionMatrix::operator()(Symbol a,Symbol b)
{
  return M(a,b);//alphabetMap(a),alphabetMap(b));
}



void SubstitutionMatrix::printOn(ostream &os) const
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
      os << alphabet.lookup(i); //alphabetMap.unmap(i));
      for(int i=1 ; i<colWidth ; ++i) os << ' ';
    }
  os << endl;

  // Print each row of the matrix
  for(Symbol i=0 ; i<N ; ++i)
    {
      os << alphabet.lookup(i)<<' ';//alphabetMap.unmap(i)) << ' ';
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



ostream &operator<<(ostream &os,const SubstitutionMatrix &M)
{
  M.printOn(os);
  return os;
}



void SubstitutionMatrix::convertToLogs()
{
  if(logSpace) throw "SubstitutionMatrix::converToLogs() : in log space";
  int n=M.getNumRows();
  for(int i=0 ; i<n ; ++i)
    for(int j=0 ; j<n ; ++j)
      M(i,j)=log(M(i,j));
  logSpace=true;
}



void SubstitutionMatrix::initGenerators() 
{
  if(logSpace) throw "SubstitutionMatrix::initGenerators() : in log space";
  int numRows=M.getNumRows();
  generators.resize(numRows);
  for(int i=0 ; i<numRows ; ++i) {
    RouletteWheel &wheel=generators[i];
    for(int j=0 ; j<numRows ; ++j)
      wheel.addSector(M(i,j));
    wheel.doneAddingSectors();
  }
}



Symbol SubstitutionMatrix::mutate(Symbol s) {
    RouletteWheel &wheel=generators[alphabetMap(s)];
    //cout<<wheel<<endl;
    Symbol r=wheel.spin();
    r=alphabetMap.unmap(r);
    return r;
}



void SubstitutionMatrix::getEqFreqs(Array1D<double> &eqFreqs)
{
  if(logSpace) throw "SubstitutionMatrix::getEqFreqs() : in log space";
  GSL::Matrix G; // eigenvectors
  GSL::Vector eigenvalues;
  GSL::Vector imagParts; // should be all zeros for practical rate matrices
  M.getEigenVectors(G,eigenvalues,imagParts);
  GSL::Matrix G_inv;
  G.invert(G_inv);
  double a=G(0,0);
  int n=M.getNumRows();
  eqFreqs.resize(n);
  for(int i=0 ; i<n ; ++i)
    eqFreqs[i]=a*G_inv(0,i);
}



void SubstitutionMatrix::setAlphabetMap(const AlphabetMap &A)
{
  alphabetMap=A;
}


