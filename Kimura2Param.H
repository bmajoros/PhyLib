/****************************************************************
 Kimura2Param.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Kimura2Param_H
#define INCL_Kimura2Param_H
#include "RateMatrix.H"

class Kimura2Param : public RateMatrix
{
    double transitionP, transversionP;
public:
    Kimura2Param(double transitionProb,double transversionProb);
    Kimura2Param(istream &is);
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

