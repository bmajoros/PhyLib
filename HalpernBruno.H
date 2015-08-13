/****************************************************************
 HalpernBruno.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_HalpernBruno_H
#define INCL_HalpernBruno_H
#include <iostream>
#include "RateMatrix.H"
#include "BOOM/Array1D.H"
#include "BOOM/GSL/ContinuousDistribution.H"
using namespace std;
using namespace BOOM;


class HalpernBruno : public RateMatrix {
  Array1D<double> eq;
  RateMatrix *Q;
public:
  HalpernBruno(RateMatrix &backgroundMutationRate,
	       const Array1D<double> &eqFreqs);
  HalpernBruno(istream &);
  virtual ~HalpernBruno();
  virtual void save(ostream &);
  virtual int numParameters() const;
  virtual double getIthParm(int i) const;
  virtual void setIthParm(int i,double parm);
  virtual RateMatrix *average(const Array1D<RateMatrix*> &) const;
  virtual void addNoise(GSL::ContinuousDistribution &);
  virtual RateMatrix *clone() const;
  virtual void init();
};

#endif

