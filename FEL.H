/****************************************************************
 FEL.H : Felsenstein's rate matrix
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_FEL_H
#define INCL_FEL_H
#include "BOOM/Array1D.H"
#include "RateMatrix.H"


class FEL : public RateMatrix
{
  BOOM::Array1D<double> eq;
  double alpha;
public:
  FEL(BOOM::Array1D<double> equilibriumFrequencies,double mutationRate);
  FEL(istream &);
  virtual RateMatrixType getType();
  virtual void save(ostream &);
  virtual int numParameters() const;
  virtual void partialDerivative(int parameter,RateMatrix &) const;
  virtual double getIthParm(int i) const;
  virtual void setIthParm(int i,double parm);
  virtual RateMatrix *average(const Array1D<RateMatrix*> &) const;
  virtual void addNoise(GSL::ContinuousDistribution &);
  virtual RateMatrix *clone() const;
  virtual void init();
};

#endif

