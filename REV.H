/****************************************************************
 REV.H : General reversible matrix (GTR)
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_REV_H
#define INCL_REV_H
#include "RateMatrix.H"
#include "BOOM/Array1D.H"


class REV /* =GTR */ : public RateMatrix
{
  BOOM::Array1D<double> eq;
  double alpha, beta, kappa, chi, omega, tau;
  void bound(double &param);
  void sumToLessThanOne(double &a,double &b,double &c);    
public:
  REV(BOOM::Array1D<double> equilibriumFrequencies,double alpha,
      double beta,double kappa,double chi,double omega,double tau);
  REV(istream &);
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
