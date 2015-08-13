#ifndef INCL_FelsensteinInterface_H
#define INCL_FelsensteinInterface_H


class FelsensteinInterface {
public:
  virtual double logLikelihood(int column)=0;
  virtual double logLikelihood()=0;
  virtual double logLikelihood(int begin,int end)=0;
  virtual double logLikelihood3(int begin,int end)=0;//3rd codon position only
  virtual double logLikelihood_bestFrame(int begin,int end)=0;
  virtual double logLikelihoodInPhase(int begin,int end,int phase)=0;
};


#endif
