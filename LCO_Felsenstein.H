/****************************************************************
 LCO_Felsenstein.H (Leaf Contexts Only)
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_LCO_Felsenstein_H
#define INCL_LCO_Felsenstein_H
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



class LCO_Felsenstein;



class LCO_FelThread : public Thread
{
    double L;
    int begin, end; // columns
    LCO_Felsenstein &alg;
public:
    LCO_FelThread(LCO_Felsenstein &,int begin,int end);
    virtual void f();
    double getLogLikelihood();
};



class LCO_Felsenstein : public FelsensteinInterface
{
public:
  LCO_Felsenstein(int order,Phylogeny &,NthOrdRateMatrix &,MultSeqAlignment&,
		    MolecularSequenceType,int numThreads=1);
  virtual double logLikelihood(int column);
  virtual double logLikelihood();
  virtual double logLikelihood(int begin,int end);
  virtual double logLikelihood_bestFrame(int begin,int end);
  virtual double logLikelihood3(int begin,int end);
  virtual double logLikelihoodInPhase(int begin,int end,int phase);
private:
  Phylogeny &phylogeny;
  NthOrdRateMatrix &Q;
  MultSeqAlignment &alignment;
  MolecularSequenceType seqType;
  const Alphabet &alphabet;
  int numNodes;
  int numAlpha;
  int largestID;
  int numThreads;
  Array1D<LCO_FelThread*> threads;
  int order;
  
  void attachPhylogenyNodes();
  double logLikelihood_threaded();
};



#endif

