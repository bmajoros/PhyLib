/****************************************************************
 GapRateMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_GapRateMatrix_H
#define INCL_GapRateMatrix_H
#include <iostream>
#include "BOOM/GSL/Matrix.H"
#include "BOOM/MolecularSequenceType.H"
#include "BOOM/GapPatternAlphaMap.H"
#include "BOOM/GSL/ContinuousDistribution.H"
#include "BOOM/GapPatternAlphabet.H"
#include "SubstitutionMatrix.H"
#include "RateMatrix.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 GapRateMatrix : a rate matrix over the alphabet {N,-,.}, denoting
 nucleotides, gaps, and unaligned positions.
 ****************************************************************/
class GapRateMatrix : public RateMatrix
{
public:
  GapRateMatrix();
  GapRateMatrix(istream &);
  virtual void save(ostream &);
  
  virtual int numParameters() const;
  virtual double getIthParm(int i) const;
  virtual void setIthParm(int i,double parm);
  virtual void partialDerivative(int parameter,GapRateMatrix &) const;
  virtual RateMatrix *average(const Array1D<RateMatrix*> &) const;
  virtual void addNoise(GSL::ContinuousDistribution &);
  virtual GapRateMatrix *clone() const;
  virtual void init();
  static GapRateMatrix *random();
protected:
  GapPatternAlphabet alphabet;
  double dashParm; // P(N -> dash)
  double dotParm;  // P(N -> dot)
  
  friend ostream &operator<<(ostream &,const GapRateMatrix &);
  void printOn(ostream &) const;
  void load(istream &);
};

ostream &operator<<(ostream &,const GapRateMatrix &);


#endif
