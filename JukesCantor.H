/****************************************************************
 JukesCantor.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_JukesCantor_H
#define INCL_JukesCantor_H
#include "RateMatrix.H"
#include "BOOM/String.H"

class JukesCantor : public RateMatrix
{
    double alpha;
public:
    JukesCantor(double alpha);
    JukesCantor(istream &);
    virtual int numParameters() const;
    virtual void partialDerivative(int parameter,RateMatrix &) const;
    virtual void save(ostream &);
    virtual RateMatrixType getType();
    virtual double getIthParm(int i) const;
    virtual void setIthParm(int i,double parm);
    virtual RateMatrix *average(const Array1D<RateMatrix*> &) const;
    virtual void addNoise(GSL::ContinuousDistribution &);
    virtual RateMatrix *clone() const;
    virtual void init();
};

#endif

