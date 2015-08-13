/****************************************************************
 FelsensteinsAlgorithm.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_FelsensteinsAlgorithm_H
#define INCL_FelsensteinsAlgorithm_H
#include <iostream>
#include "BOOM/MultiAlignment.H"
#include "BOOM/Alphabet.H"
#include "BOOM/MolecularSequenceType.H"
#include "BOOM/Vector.H"
#include "BOOM/Thread.H"
#include "Phylogeny.H"
#include "RateMatrix.H"
#include "FelsensteinInterface.H"
using namespace std;
using namespace BOOM;



class FelsensteinsAlgorithm;



class FelThread : public Thread
{
    double L;
    int begin, end; // columns
    FelsensteinsAlgorithm &alg;
public:
    FelThread(FelsensteinsAlgorithm &,int begin,int end);
    virtual void f();
    double getLogLikelihood();
};



class FelsensteinsAlgorithm : public FelsensteinInterface
{
public:
  FelsensteinsAlgorithm(Phylogeny &,RateMatrix &,MultiAlignment &,
			MolecularSequenceType,int numThreads=1);
  virtual double logLikelihood(int column);
  virtual double logLikelihood();
  virtual double logLikelihood(int begin,int end);
private:
    Phylogeny &phylogeny;
    RateMatrix &Q;
    MultiAlignment &alignment;
    MolecularSequenceType seqType;
    Alphabet &alphabet;
    int numNodes;
    int numAlpha;
    int largestID;
    int numThreads;
    Array1D<FelThread*> threads;

    void attachPhylogenyNodes();
    double logLikelihood_threaded();
};



#endif

