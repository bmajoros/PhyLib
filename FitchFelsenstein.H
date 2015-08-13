/****************************************************************
 FitchFelsenstein.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_FitchFelsenstein_H
#define INCL_FitchFelsenstein_H
#include <iostream>
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Alphabet.H"
#include "BOOM/MolecularSequenceType.H"
#include "BOOM/Vector.H"
#include "BOOM/Thread.H"
#include "Phylogeny.H"
#include "NthOrdRateMatrix.H"
#include "FelsensteinInterface.H"
using namespace std;
using namespace BOOM;



class FitchFelsenstein;



class FitchThread : public Thread
{
    double L;
    int begin, end; // columns
    FitchFelsenstein &alg;
public:
    FitchThread(FitchFelsenstein &,int begin,int end);
    virtual void f();
    double getLogLikelihood();
};



class FitchFelsenstein : public FelsensteinInterface
{
public:
  FitchFelsenstein(int order,Phylogeny &,NthOrdRateMatrix &,MultSeqAlignment&,
		    MolecularSequenceType,int numThreads=1);
  virtual double logLikelihood(int column);
  virtual double logLikelihood();
  virtual double logLikelihood(int begin,int end);
  virtual double logLikelihood_bestFrame(int begin,int end);
  virtual double logLikelihood3(int begin,int end);//3rd codon position only
  virtual double logLikelihoodInPhase(int begin,int end,int phase);
private:
  Phylogeny &phylogeny;
  NthOrdRateMatrix &Q;
  MultSeqAlignment &alignment;
  MolecularSequenceType seqType;
  Alphabet &alphabet;
  int numNodes;
  int numAlpha;
  int largestID;
  int numThreads;
  Array1D<FitchThread*> threads;
  int order;
  
  void attachPhylogenyNodes();
  double logLikelihood_threaded();
};



#endif

