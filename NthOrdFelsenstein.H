/****************************************************************
 NthOrdFelsenstein.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_NthOrdFelsenstein_H
#define INCL_NthOrdFelsenstein_H
#include <iostream>
#include "BOOM/Alphabet.H"
#include "ContextType.H"
#include "Phylogeny.H"
#include "NthOrdRateMatrix.H"
#include "BOOM/MolecularSequenceType.H"
#include "BOOM/MultSeqAlignment.H"
#include "FelsensteinInterface.H"
#include "AlignmentNmerTable.H"
using namespace std;
using namespace BOOM;


class NthOrdFelsenstein : public FelsensteinInterface
{
public:
  NthOrdFelsenstein(ContextType,Phylogeny &,NthOrdRateMatrix *,
		    MultSeqAlignment &,MolecularSequenceType,
		    AlignmentNmerTable *,int numThreads=1);
  virtual double logLikelihood(int column);
  virtual double logLikelihood();
  virtual double logLikelihood(int begin,int end);
  virtual double logLikelihood3(int begin,int end);//3rd codon positions only
  virtual double logLikelihoodInPhase(int begin,int end,int phase);
  virtual double logLikelihood_bestFrame(int begin,int end) 
    {throw "NthOrdFelsenstein::logLikelihood_bestFrame()";}
private:
  ContextType contextType;
  Phylogeny &phylogeny;
  NthOrdRateMatrix *Q;
  MultSeqAlignment &A;
  MolecularSequenceType seqType;
  int numThreads;
  Alphabet &alphabet;
  int order;
  bool threePeriodic;
  AlignmentNmerTable *nmerTable;
};

#endif

