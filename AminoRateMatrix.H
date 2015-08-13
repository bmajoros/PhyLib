/****************************************************************
 AminoRateMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_AminoRateMatrix_H
#define INCL_AminoRateMatrix_H
#include <iostream>
#include "RateMatrix.H"
using namespace std;
using namespace BOOM;


class AminoRateMatrix : public RateMatrix
{
public:
    AminoRateMatrix();
    virtual void save(ostream &);
    virtual RateMatrixType getType();
    virtual int numParameters() const;
    virtual double getIthParm(int i) const;
    virtual void setIthParm(int i,double parm);
    virtual void partialDerivative(int parameter,RateMatrix &) const;
    virtual RateMatrix *average(const Array1D<RateMatrix*> &) const;
    virtual void addNoise(GSL::ContinuousDistribution &);
    virtual RateMatrix *clone() const;
};

#endif

