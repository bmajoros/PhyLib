/****************************************************************
 MPC.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_MPC_H
#define INCL_MPC_H
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



class MPC;



class FitchThread : public Thread
{
    double L;
    int begin, end; // columns
    MPC &alg;
public:
    FitchThread(MPC &,int begin,int end);
    virtual void f();
    double getLogLikelihood();
};



class MPC
{
public:
  MPC(int order,Phylogeny &,NthOrdRateMatrix &,MultSeqAlignment&,
		    MolecularSequenceType,int numThreads=1);
  virtual double logLikelihood(int column);
  virtual double logLikelihood();
  virtual double logLikelihood(int begin,int end);
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

