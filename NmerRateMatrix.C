/****************************************************************
 NmerRateMatrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include <math.h>
#include "NmerRateMatrix.H"
#include "NmerSubstMatrix.H"
#include "BOOM/Random.H"
#include "BOOM/Constants.H"
#include "BOOM/PureDnaAlphabet.H"
using namespace std;
using namespace BOOM;


NmerRateMatrix::NmerRateMatrix(int n)
  : NthOrdRateMatrix(n-1,DNA,MT_NMER,true),
    eqSingle(4), raw(false)
{
  // ctor

  // Initialize the parameter vector
  int numParms=numParameters();
  parms.resize(numParms);
  parms.setAllTo(0.01); // in case the user forgets to set them later

  // Initialize the equilibrium frequencies
  numNmers=pow((float)alphabetMap.getRangeSize(),(float)n);
  eq.resize(numNmers);
  float p=1.0/numNmers;
  eq.setAllTo(p);
  eqSingle.setAllTo(0.25);

  // Initialize the mapping from nmer pairs to parameter indices
  initParmIndices();

  // Allocate memory for the matrix
  matrix.resize(numNmers,numNmers);

  // Initialize any lower order models if applicable
  if(n>1) {
    delete lowerOrderModel;
    lowerOrderModel=new NmerRateMatrix(n-1);
  }
}



NmerRateMatrix::~NmerRateMatrix()
{
  // dtor
}



const NmerRateMatrix &NmerRateMatrix::operator=(const NmerRateMatrix &other)
{
  parms=other.parms;
  eq=other.eq;
  eqSingle=other.eqSingle;
  matrix=other.matrix;
  raw=other.raw;

  if(lowerOrderModel)
    *lowerOrderModel=*other.lowerOrderModel;

  // Base class has some attributes too, but I don't care about those...
  
  return *this;
}



void NmerRateMatrix::eqFreqsFromMarginals(Array1D<double> singleNucEqFreqs)
{
  eqSingle=singleNucEqFreqs;

  const int n=order+1;
  Sequence nmer;
  for(int i=0 ; i<numNmers ; ++i) {
    nmer.fromInt(i,n,alphabetMap);
    double P=1.0;
    for(int pos=0 ; pos<n ; ++pos) {
      P*=singleNucEqFreqs[nmer[pos]]; // ### should be in log space??
      //cout<<singleNucEqFreqs[nmer[pos]]<<" * ";
    }
    //cout<<" = "<<P<<endl;
    eq[i]=P;
  }

  if(lowerOrderModel)
    static_cast<NmerRateMatrix*>(lowerOrderModel)->
      eqFreqsFromMarginals(singleNucEqFreqs);
}



void NmerRateMatrix::randomize(double noiseSigma,
			       Array1D<double> singleNucEqFreqs)
{
  eqFreqsFromMarginals(singleNucEqFreqs);
  int n=parms.size();
  for(int i=0 ; i<n ; ++i)
    parms[i]=RandomFloat(0.01,0.10);

  if(lowerOrderModel) 
    static_cast<NmerRateMatrix*>(lowerOrderModel)->
      randomize(noiseSigma,singleNucEqFreqs);
}



void NmerRateMatrix::setNmerEqFreq(const Sequence &nmer,double nmerFreq)
{
  int index=nmer.asInt(alphabetMap);
  eq[index]=nmerFreq;
}



int NmerRateMatrix::numParameters() const
{
  const int n=order+1;
  return 6*n*pow(float(alphabetMap.getRangeSize()),float(n-1));
}



double NmerRateMatrix::getIthParm(int i) const
{
  return parms[i];
}



void NmerRateMatrix::setIthParm(int i,double parm)
{
  parms[i]=parm;
}



int NmerRateMatrix::countDifferences(const Sequence &seq1,
				     const Sequence &seq2)
{
  int length=seq1.getLength(), diffs=0;
  for(int i=0 ; i<length ; ++i)
    if(seq1[i]!=seq2[i])
      ++diffs;
  return diffs;
}



void NmerRateMatrix::initParmIndices()
{
  const int n=order+1;
  Sequence from, to;
  int nextParm=0;
  for(int i=0 ; i<numNmers ; ++i) {
    from.fromInt(i,n,alphabetMap);
    //cout<<i<<"=";from.printOn(cout,*alphabetMap.getRange());cout<<endl;
    for(int j=i+1 ; j<numNmers ; ++j) {
      to.fromInt(j,n,alphabetMap);
      int numDifferences=countDifferences(from,to);
      if(numDifferences!=1) continue;
      Sequence from_to=from+to; 
      Sequence to_from=to+from; 
      parmIndices[from_to]=parmIndices[to_from]=nextParm;
      //cout<<nextParm<<"=";from_to.printOn(cout,*alphabetMap.getRange());cout<<endl;
      ++nextParm;
    }
  }
}



void NmerRateMatrix::init()
{
  if(raw) return;
  matrix.setAllTo(0.0);
  const int n=order+1;
  Sequence from, to;
  int nextParm=0;
  for(int i=0 ; i<numNmers ; ++i) {
    from.fromInt(i,n,alphabetMap);
    for(int j=i+1 ; j<numNmers ; ++j) {
      to.fromInt(j,n,alphabetMap);
      int numDifferences=countDifferences(from,to);
      if(numDifferences!=1) continue;
      double p=parms[nextParm];
      matrix(i,j)=p*eq[j];
      matrix(j,i)=p*eq[i];
      ++nextParm;
    }
  }

  // Install diagonal to bring row sums to zero
  for(int i=0 ; i<numNmers ; ++i) {
    double sum=0.0;
    for(int j=0 ; j<numNmers ; ++j)
      sum+=matrix(i,j);
    matrix(i,i)=-sum;
  }

  // Scale entire matrix to a constant substitution frequency
  // (Yang & Nielsen 2000, MBE vol. 17, p33, Eq. #3)
  /* COMMENTED OUT -- FIXING ONE BRANCH LENGTH INSTEAD
  double sum=0;
  for(int i=0 ; i<numNmers ; ++i) {
    double innerSum=0;
    for(int j=0 ; j<numNmers ; ++j)
      if(j!=i) innerSum+=matrix(i,j);
    sum+=eq[i]*innerSum;
  }
  double scalingFactor=1/sum;
  for(int i=0 ; i<numNmers ; ++i)
    for(int j=0 ; j<numNmers ; ++j)
      matrix(i,j)*=scalingFactor;
  */
}



int NmerRateMatrix::getNumRows() const 
{
  throw "NmerRateMatrix::getNumRows()";
}



void NmerRateMatrix::averageLowerOrderModel(bool useDualContexts)
{
  if(!lowerOrderModel) return;
  NmerRateMatrix *M=static_cast<NmerRateMatrix*>(lowerOrderModel);
  M->averageFromHigherOrderModel(this);
}



void NmerRateMatrix::averageFromHigherOrderModel(NmerRateMatrix *H)
{
  Sequence from, to;
  int n=order+1;
  for(int i=0 ; i<numNmers ; ++i) {
    from.fromInt(i,n,alphabetMap);
    for(int j=i+1 ; j<numNmers ; ++j) {
      to.fromInt(j,n,alphabetMap);
      int numDifferences=countDifferences(from,to);
      if(numDifferences!=1) continue;
      Sequence key=from+to;
      int parmIndex=parmIndices[key];
      double &myParm=parms[parmIndex];
      double sum=0.0;
      int sampleSize=0;
      for(Symbol x=0 ; x<4 ; ++x) {//###
	Sequence hiFrom=from+x;
	int fromCode=hiFrom.asInt(alphabetMap);
	Sequence hiTo=to+x;
	int toCode=hiTo.asInt(alphabetMap);
	//if(toCode<fromCode) continue;
	key=hiFrom+hiTo;
	if(!H->parmIndices.isDefined(key)) continue;
	parmIndex=H->parmIndices[key];
	sum+=eqSingle[x]*H->parms[parmIndex];//eqSingle[x]*eqSingle[y]*H->parms[parmIndex];
	++sampleSize;
      }
      myParm=sum;//sum/sampleSize; //sum;
    }
  }
  averageLowerOrderModel(false);
}



int NmerRateMatrix::lookupIndex(const Sequence &parentContext,
				const Sequence &childContext)
{
  throw "NmerRateMatrix::lookupIndex";
}



int NmerRateMatrix::getNumMatrices() const
{
  return 1;
}



NthOrdRateMatrix *NmerRateMatrix::nextHigherOrder()
{
  // Misc initialization
  const int n=order+1; // my n, not M's!
  const int Nm=n+1;    // M's n
  NmerRateMatrix *M=new NmerRateMatrix(Nm);
  Sequence from, to, fromSuffix, toSuffix, fromPrefix, toPrefix;
  int numNmersInM=M->numNmers;
  M->eqFreqsFromMarginals(eqSingle);
  *(M->lowerOrderModel)=*this;
  M->parms.setAllTo(-1);

  // Compute averages of parameter values for each nmer in the lower-order
  // model, for use in initializing higher-order parameters having no direct
  // analog in the lower-order model
  Array1D<double> aveParms(numNmers);
  for(int i=0 ; i<numNmers ; ++i) {
    from.fromInt(i,n,alphabetMap);
    double sum=0;
    int sampleSize=0;
    for(int j=0 ; j<numNmers ; ++j) {
      to.fromInt(j,n,alphabetMap);
      Sequence fromTo=from+to;
      if(parmIndices.isDefined(fromTo)) {
	sum+=parms[parmIndices[fromTo]];
	++sampleSize;
      }
    }
    aveParms[i]=1-sum; //  sum/sampleSize;
  }

  NmerRateMatrix *zerothM=this;
  if(order>0) zerothM=(NmerRateMatrix*)getLowerOrderModel(0);
  //if(!zerothM) throw "zerothM is null";
  int numAlpha=alphabetMap.getDomain()->size();
  Array1D<double> zerothParm(numAlpha);
  zerothParm.setAllTo(NEGATIVE_INFINITY);
  Sequence seq;
  seq.resize(2);
  for(Symbol s=0 ; s<numAlpha ; ++s) {
    if(alphabetMap(s)==INVALID_SYMBOL) continue;
    zerothParm[s]=1;
    for(Symbol z=0 ; z<numAlpha ; ++z) {
      if(z==s || alphabetMap(z)==INVALID_SYMBOL) continue;
      seq[0]=s; seq[1]=z;
      zerothParm[s]-=zerothM->parms[zerothM->parmIndices[seq]];
    }
  }

  // Now iterate over all nmer pairs in the higher-order model and
  // initialize parameters from lower-order analogs (or from aveParms[])
  for(int i=0 ; i<numNmersInM ; ++i) {
    from.fromInt(i,Nm,alphabetMap);
    from.getSubsequence(1,n,fromSuffix);
    int suffixInt=fromSuffix.asInt(alphabetMap);
    for(int j=i+1 ; j<numNmersInM ; ++j) {
      to.fromInt(j,Nm,alphabetMap);
      if(countDifferences(from,to)!=1) continue;
      Sequence nmerPair=from+to;
      int higherParmIndex=M->parmIndices[nmerPair];
      if(from[0]==to[0]) {
	// The difference is in the suffix, so just copy the parameter
	// from the lower-order model
	to.getSubsequence(1,n,toSuffix);
	nmerPair=fromSuffix+toSuffix;
	int lowerParmIndex=parmIndices[nmerPair];
	M->parms[higherParmIndex]=
	  parms[lowerParmIndex] *
	  zerothParm[from[0]];
      }
      /*else {
	// The difference is in the first position, so copy the parameter
	// from the lower-order model for the prefixes (not suffixes)
	from.getSubsequence(0,n,fromPrefix);
	to.getSubsequence(0,n,toPrefix);
	nmerPair=fromPrefix+toPrefix;
	int lowerParmIndex=parmIndices[nmerPair];
	M->parms[higherParmIndex]=parms[lowerParmIndex]
	  ;//*zerothParm[from[Nm-1]];
	  }*/
      else {
	// The difference is in the leftmost base, so just initialize
	// the new parameter to an average of the parameters for the
	// row of the suffix in the lower-order model
	seq[0]=from[0]; seq[1]=to[0];
	M->parms[higherParmIndex]=
	  aveParms[suffixInt] *
	  zerothM->parms[zerothM->parmIndices[seq]];
      }
    }
  }

  return M;
}



void NmerRateMatrix::save(const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  save(os);
}



void NmerRateMatrix::save(ostream &os)
{
  os<<NOMT_NMER<<endl; // NthOrdMatrixType
  os<<order<<endl;
  os<<parms<<endl;
  os<<eq<<endl;
  os<<eqSingle<<endl;
  
  if(lowerOrderModel)
    lowerOrderModel->save(os);
}



unsigned convertNmerCode(unsigned rawCode,int seqLength,const Alphabet &alphabet)
{
  Sequence S;
  S.resize(seqLength);
  unsigned base=alphabet.size();
  for(int i=seqLength-1 ; i>=0 ; --i) {
    unsigned digit=rawCode%base;
    Symbol s=(int)digit;
    S[i]=s;
    rawCode/=base;
  }	
  return S.asInt(alphabet,0,S.getLength());
}



NthOrdRateMatrix *NmerRateMatrix::loadRaw(istream &is)
{
  // Reads from a raw matrix file

  int order;
  is>>order;
  NmerRateMatrix *Q=new NmerRateMatrix(order+1);
  is>>Q->eq;

  //is>>Q->matrix;
  Alphabet alphabet=PureDnaAlphabet::global();
  GSL::Matrix M;
  is>>M;
  unsigned n=M.getNumRows();
  for(unsigned i=0 ; i<n ; ++i) {
    unsigned newI=convertNmerCode(i,order+1,alphabet);
    for(unsigned j=0 ; j<n ; ++j) {
      unsigned newJ=convertNmerCode(j,order+1,alphabet);
      Q->matrix(newI,newJ)=M(i,j);
    }
  }


  Q->raw=true;

  // ### still need to infer eqSingle
  
  if(order>0) {
    NmerRateMatrix *m=new NmerRateMatrix(order);
    //m->averageFromHigherOrderModel(Q);
    Q->lowerOrderModel=m;
  }
  return Q;
}



NthOrdRateMatrix *NmerRateMatrix::load(istream &is)
{
  // Precondition: the type indicator (NthOrdMatrixType) has already 
  // beeen read, and is known to be NOMT_NMER

  int order;
  is>>order;
  NmerRateMatrix *Q=new NmerRateMatrix(order+1);
  is>>Q->parms;
  is>>Q->eq;
  is>>Q->eqSingle;

  if(order>0)
    Q->lowerOrderModel=NthOrdRateMatrix::load(is);

  return Q;
}



RateMatrix &NmerRateMatrix::lookup(const Sequence &context,int phase,
				   int begin,int len)
{
  throw "NmerRateMatrix::lookup()";
}



RateMatrix &NmerRateMatrix::lookup(const Sequence &parentContext,
				   const Sequence &childContext,
				   int phase)
{
  throw "NmerRateMatrix::lookup()";
}



RateMatrix &NmerRateMatrix::getIthMatrix(int i,int phase)
{
  throw "NmerRateMatrix::lookup()";
}



void NmerRateMatrix::instantiate(double branchLength,NthOrdSubstMatrix &pt)
{
  init();

  NmerSubstMatrix *Pt=static_cast<NmerSubstMatrix*>(&pt);
  rateMatrixToSubstMatrix(matrix,branchLength,Pt->peek());

  //Pt->averageLowerOrderModel(eqSingle);
  Pt->convertToLogs();

  //NmerRateMatrix *Q=static_cast<NmerRateMatrix*>(lowerOrderModel);
  //if(Q) Q->instantiate(branchLength,*Pt->getLowerOrderModel());
}



NthOrdSubstMatrix *NmerRateMatrix::instantiate(double branchLength)
{
  NthOrdSubstMatrix *Pt=new NmerSubstMatrix(order,alphabetMap);
  instantiate(branchLength,*Pt);
  return Pt;
}



NthOrdRateMatrix *NmerRateMatrix::ancestralContextsOnly()
{
  throw "NmerRateMatrix::ancestralContextsOnly()";
}



NthOrdRateMatrix *NmerRateMatrix::clone()
{
  return new NthOrdRateMatrix(*this);
}



bool NmerRateMatrix::isPeriodic()
{
  return false;
}



void NmerRateMatrix::printOn(ostream &os) const
{
  os<<matrix<<endl<<endl;
  parms.printOn(os);
  os<<endl<<endl;
  eq.printOn(os);
  os<<endl<<endl;
}



double NmerRateMatrix::getNmerEqFreq(int nmerCode)
{
  return eq[nmerCode];
}


