/****************************************************************
 NmerRateMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_NmerRateMatrix_H
#define INCL_NmerRateMatrix_H
#include "NthOrdRateMatrix.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/Map.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 NmerRateMatrix : a rate matrix describing instantaneous rates of
 substitution from one nmer to another nmer. In this class we
 assume that at most one position in the nmer changes; any pair
 of nmers differing in more than one position will have a rate of
 zero.
 ****************************************************************/
class NmerRateMatrix : public NthOrdRateMatrix
{
public:
  NmerRateMatrix(int N); // note! N=order+1
  ~NmerRateMatrix();
  const NmerRateMatrix &operator=(const NmerRateMatrix &);
  static NthOrdRateMatrix *load(istream &);
  static NthOrdRateMatrix *loadRaw(istream &);
  void setNmerEqFreq(const Sequence &nmer,double nmerfreq);
  double getNmerEqFreq(int nmerCode);
  void randomize(double noiseSigma,Array1D<double> singleNucEqFreqs);
  virtual void init();
  virtual int getNumRows() const;
  virtual void averageLowerOrderModel(bool useDualContexts=false);
  virtual int lookupIndex(const Sequence &parentContext,
			  const Sequence &childContext);
  virtual void instantiate(double branchLength,NthOrdSubstMatrix &);
  virtual int getNumMatrices() const;
  virtual int numParameters() const;
  virtual NthOrdRateMatrix *nextHigherOrder();
  virtual void save(const BOOM::String &filename);
  virtual void save(ostream &);
  virtual RateMatrix &lookup(const Sequence &context,int phase=-1,int begin=0,
			     int len=-1);
  virtual RateMatrix &lookup(const Sequence &parentContext,
			     const Sequence &childContext,
			     int phase=-1);
  virtual RateMatrix &getIthMatrix(int i,int phase=-1);
  virtual NthOrdSubstMatrix *instantiate(double branchLength);
  virtual NthOrdRateMatrix *ancestralContextsOnly();
  virtual NthOrdRateMatrix *clone();
  virtual bool isPeriodic();
  virtual double getIthParm(int i) const;
  virtual void setIthParm(int i,double parm);
  virtual void printOn(ostream &) const;
  
protected:
  int numNmers;
  Array1D<double> parms;
  Array1D<double> eq; // nmer equilibrium frequencies
  Array1D<double> eqSingle; // single-nucleotide equilibrium freqs
  GSL::Matrix matrix;
  BOOM::Map<Sequence,int> parmIndices; // concatenate sequences to get the key
  bool raw;

  void initParmIndices();
  int countDifferences(const Sequence &,const Sequence &);
  void eqFreqsFromMarginals(Array1D<double> singleNucEqFreqs);
  void averageFromHigherOrderModel(NmerRateMatrix *higherOrderModel);
};

#endif

