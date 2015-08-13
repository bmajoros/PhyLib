/****************************************************************
 NthOrdRateMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_NthOrdRateMatrix_H
#define INCL_NthOrdRateMatrix_H
#include <iostream>
#include "RateMatrix.H"
#include "BOOM/Array1D.H"
#include "BOOM/Sequence.H"
#include "NthOrdSubstMatrix.H"
using namespace std;
using namespace BOOM;



enum NthOrdMatrixType {
  NOMT_COND, // NthOrdRateMatrix ("conditional")
  NOMT_NMER, // NmerRateMatrix ("joint", but w/only 1 mutation per generation)
  NOMT_NMER_RAW // same as NMER, but read from a raw file (no parms)
};
ostream &operator<<(ostream &,NthOrdMatrixType);
istream &operator>>(istream &,NthOrdMatrixType &);



/****************************************************************
 This class represents higher-order rate matrices, which are 
 conditional on some number of preceding residues in the ancestor 
 species (see descendent class NmerRateMatrix for the use of dual 
 contexts).
 ****************************************************************/
class NthOrdRateMatrix
{
public:
  // Constructors/destructors: ---------------------------------------------
  NthOrdRateMatrix(int order,MolecularSequenceType,RateMatrixType,bool dual);
  virtual ~NthOrdRateMatrix();

  // Static methods: -------------------------------------------------------
  static NthOrdRateMatrix *load(istream &);
  static NthOrdRateMatrix *load(const String &filename);
  static NthOrdRateMatrix *random(bool dual,
				  double noiseSigma,
				  MolecularSequenceType,int order,
				  RateMatrixType,
				  BOOM::Array1D<double> 
				  equilibriumFreqs=omitted);

  // Instance methods: -----------------------------------------------------
  RateMatrixType getType();
  MolecularSequenceType getSeqType() {return seqType;}
  bool isDual(); // uses dual contexts? (parent & child)
  int getOrder() const;
  const Alphabet &getAlphabet() const;
  RateMatrixType getMatrixType();
  AlphabetMap &getAlphabetMap();
  const NthOrdRateMatrix &operator=(const NthOrdRateMatrix &);
  void replaceIthMatrix(int i,RateMatrix *);

  // Virtual methods: ------------------------------------------------------
  virtual void init();
  void rescale(); // called by init()
  virtual int getNumRows() const;
  virtual void averageLowerOrderModel(bool useDualContexts=false);
  virtual int lookupIndex(const Sequence &parentContext,
			  const Sequence &childContext);
  virtual void instantiate(double branchLength,NthOrdSubstMatrix &);
  virtual int getNumMatrices() const;
  virtual int numParameters() const;
  virtual Array1D<unsigned> *getParmToContextMap();
  virtual NthOrdRateMatrix *nextHigherOrder();
  virtual void save(const BOOM::String &filename);
  virtual void save(ostream &);
  virtual RateMatrix &lookup(unsigned contextCode,int contextLength);
  virtual RateMatrix &lookup(const Sequence &context,int phase=-1,int begin=0,
			     int len=-1);
  virtual RateMatrix &lookup(const Sequence &parentContext,
			     const Sequence &childContext,
			     int phase=-1);
  virtual RateMatrix &getIthMatrix(int i,int phase=-1);
  virtual NthOrdSubstMatrix *instantiate(double branchLength);
  virtual NthOrdRateMatrix *getLowerOrderModel(int order=-1,int phase=-1);
  virtual NthOrdRateMatrix *ancestralContextsOnly();
  virtual NthOrdRateMatrix *clone();
  virtual bool isPeriodic();
  virtual void printOn(ostream &) const;

protected:
  // Protected data members: -----------------------------------------------
  RateMatrixType matrixType;
  static BOOM::Array1D<double> omitted;
  Alphabet &alphabet;
  AlphabetMap alphabetMap;
  MolecularSequenceType seqType;
  int order; // order of parent + order of child!  (when dual)
  NthOrdRateMatrix *lowerOrderModel;
  Array1D<RateMatrix*> matrices;
  bool dual; // whether we use dual contexts (parent & child)
  BOOM::Array1D<double> eqFreqs;
  Array1D<unsigned> *parmToContextMap;

  // Protected methods: ----------------------------------------------------
  void randomize(double noiseSigma,BOOM::Array1D<double> equilibriumFreqs);
  void randomizeDual(double noiseSigma,
		     BOOM::Array1D<double> equilibriumFreqs);
  //friend ostream &operator<<(ostream &,const NthOrdRateMatrix &);
  //virtual void printOn(ostream &) const;
  void averageDual(); // ensures lower orders are consistent with higher
  void deleteOddOrders(); // for use in dual-context model
  void ancestralContextsOnly(NthOrdRateMatrix &);
  void clone(NthOrdRateMatrix &);
  NthOrdRateMatrix *nextHigherOrder_dual();
  double computeContextPrior(const Sequence &context);
};

ostream &operator<<(ostream &,const NthOrdRateMatrix &);


#endif
