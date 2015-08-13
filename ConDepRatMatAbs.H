/****************************************************************
 ConDepRatMatAbs.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ConDepRatMatAbs_H
#define INCL_ConDepRatMatAbs_H
#include <iostream>
#include "RateMatrix.H"
#include "BOOM/Array1D.H"
#include "BOOM/Sequence.H"
#include "NthOrdSubstMatrix.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 This class represents higher-order rate matrices, which are 
 conditional on some number of preceding residues in the ancestor 
 species.  THIS CLASS IS ABSTRACT.
 ****************************************************************/
class ConDepRatMatAbs {
public:
  // Implemented here:
  static ConDepRatMatAbs *load(const BOOM::String &);
  static ConDepRatMatAbs *load(istream &);

  // Pure virtual:
  virtual ConDepRatMatAbs *nextHigherOrder()=0;
  virtual void save(const BOOM::String &filename)=0;
  virtual void save(ostream &)=0;
  virtual RateMatrixType getType()=0;
  virtual bool isDual()=0; // uses dual contexts? (parent & child)
  virtual int getOrder() const=0;
  virtual RateMatrix &lookup(const Sequence &context,int begin=0,int len=-1)=0;
  virtual RateMatrix &lookup(const Sequence &parentContext,
			     const Sequence &childContext)=0;
  virtual int lookupIndex(const Sequence &parentContext,
			  const Sequence &childContext)=0;
  virtual int getNumMatrices() const=0;
  virtual RateMatrix &getIthMatrix(int)=0;
  virtual ConDepSubMatAbs *instantiate(double branchLength)=0;
  virtual void instantiate(double branchLength, &);
    int getNumRows() const;
    const Alphabet &getAlphabet() const;
    int numParameters() const;
    RateMatrixType getMatrixType();
    ConDepRatMatAbs *getLowerOrderModel(int order=-1);
    void averageLowerOrderModel(bool useDualContexts=false);
    AlphabetMap &getAlphabetMap() {return alphabetMap;}
    ConDepRatMatAbs *ancestralContextsOnly();
    ConDepRatMatAbs *clone();
    virtual bool isPeriodic() {return false;}
protected:
    RateMatrixType matrixType;
    static BOOM::Array1D<double> omitted;
    Alphabet &alphabet;
    AlphabetMap alphabetMap;
    MolecularSequenceType seqType;
    int order; // order of parent + order of child!  (when dual)
    ConDepRatMatAbs *lowerOrderModel;
    Array1D<RateMatrix*> matrices;
    bool dual; // whether we use dual contexts (parent & child)
    
    //void randomize(BOOM::Array1D<double>);
    void randomize(double noiseSigma,BOOM::Array1D<double> equilibriumFreqs);
    void randomizeDual(double noiseSigma,
                       BOOM::Array1D<double> equilibriumFreqs);
    friend ostream &operator<<(ostream &,const ConDepRatMatAbs &);
    void printOn(ostream &) const;
    void averageDual(); // ensures lower orders are consistent with higher
    void deleteOddOrders(); // for use in dual-context model
    void ancestralContextsOnly(ConDepRatMatAbs &);
    void clone(ConDepRatMatAbs &);
    ConDepRatMatAbs *nextHigherOrder_dual();
};

ostream &operator<<(ostream &,const SubstitutionMatrix &);


#endif
