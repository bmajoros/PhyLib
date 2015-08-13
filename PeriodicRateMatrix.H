/****************************************************************
 PeriodicRateMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_PeriodicRateMatrix_H
#define INCL_PeriodicRateMatrix_H
#include <iostream>
#include "NthOrdRateMatrix.H"
#include "BOOM/Array1D.H"
using namespace std;
using namespace BOOM;



class PeriodicRateMatrix : public NthOrdRateMatrix
{
public:
  PeriodicRateMatrix(int order,MolecularSequenceType,RateMatrixType,
		     bool dual);
  virtual ~PeriodicRateMatrix();
  static NthOrdRateMatrix *random(bool dual,
				  double noiseSigma,
				  MolecularSequenceType,int order,
				  RateMatrixType,
				  BOOM::Array1D<double> 
				  equilibriumFreqs=NthOrdRateMatrix::omitted);

  virtual NthOrdRateMatrix *nextHigherOrder(); // caller should delete me then
  virtual void save(const BOOM::String &filename);
  virtual void save(ostream &);
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
  void load_noIndicator(istream &);
private:
  Array1D<NthOrdRateMatrix*> matrices; // array of 3
};

#endif

