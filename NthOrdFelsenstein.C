/****************************************************************
 NthOrdFelsenstein.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "NthOrdFelsenstein.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/AminoAlphabet.H"
#include "TRCO_Felsenstein.H"
#include "RCO_Felsenstein.H"
#include "LCO_Felsenstein.H"
#include "HOG_Felsenstein.H"
#include "ACO_Felsenstein.H"
#include "FitchFelsenstein.H"
#include "FitchParsimony.H"
using namespace std;
using namespace BOOM;


NthOrdFelsenstein::NthOrdFelsenstein(ContextType contextType,
				     Phylogeny &phylogeny,
				     NthOrdRateMatrix *Q,
				     MultSeqAlignment &A,
				     MolecularSequenceType seqType,
				     AlignmentNmerTable *nmerTable,
				     int numThreads) 
  : contextType(contextType), phylogeny(phylogeny), Q(Q), A(A),
    seqType(seqType), numThreads(numThreads), order(Q->getOrder()),
    alphabet(seqType==DNA ? 
	     //(Alphabet&) DnaAlphabet::global() : 
	     (Alphabet&) DnaDashDotAlphabet::global() :
	     (Alphabet&) AminoAlphabet::global()),
    threePeriodic(Q->isPeriodic()),
    nmerTable(nmerTable)
{
  // ctor

  //if(contextType==CT_ACO) contextType=CT_RCO;

  switch(contextType) 
    {
    case CT_MP: // maximum parsimony
      {
	if(Q->isDual()) {
	  this->Q=Q->ancestralContextsOnly();
	  order/=2;
	}
	//FitchParsimony fitch(phylogeny,alphabet,A,A.getGapSymbols());
	//fitch.run();
	//phylogeny.attachAlignment(A);
      }
      break;
    case CT_TRCO: // transitive root contexts only
      {
	if(Q->isDual()) {
	  this->Q=Q->ancestralContextsOnly();
	  order/=2;
	}
	//phylogeny.attachAlignment(A);
      }
      break;
    case CT_RCO: // root contexts only
      {
	if(Q->isDual()) {
	  this->Q=Q->ancestralContextsOnly();
	  order/=2;
	}
	//phylogeny.attachAlignment(A);
      }
      break;
    case CT_LCO: // leaf contexts only
      {
	if(Q->isDual()) {
	  this->Q=Q->ancestralContextsOnly();
	  order/=2;
	}
	//phylogeny.attachAlignment(A);
      }
      break;
    case CT_ACO: // ancestral contexts only
      {
	if(Q->isDual()) {
	  this->Q=Q->ancestralContextsOnly();
	  order/=2;
	}
	//phylogeny.attachAlignment(A);
      }
      break;
    case CT_HOG: // the "whole hog" (dual contexts)
      {
	if(!Q->isDual()) throw "HOG model requires dual contexts";
	//phylogeny.attachAlignment(A);
      }
      break;
    default:
      throw "Context type not supported";
    }
}



double NthOrdFelsenstein::logLikelihood(int column) 
{
  return logLikelihood(column,column+1);
}



double NthOrdFelsenstein::logLikelihood() 
{
  return logLikelihood(0,A.getLength());
}



double NthOrdFelsenstein::logLikelihood(int begin,int end) 
{
  double LL=0;
  switch(contextType) 
    {
    case CT_MP: 
      {
	FitchFelsenstein fel(order,phylogeny,*Q,A,seqType,numThreads);
	LL=fel.logLikelihood(begin,end);
      }
      break;
    case CT_TRCO: 
      {
	TRCO_Felsenstein fel(order,phylogeny,*Q,A,seqType,numThreads);
	LL=fel.logLikelihood(begin,end);
      }
      break;
    case CT_RCO: 
      {
	RCO_Felsenstein fel(order,phylogeny,*Q,A,seqType,numThreads);
	LL=fel.logLikelihood(begin,end);
      }
      break;
    case CT_LCO: 
      {
	LCO_Felsenstein fel(order,phylogeny,*Q,A,seqType,numThreads);
	LL=fel.logLikelihood(begin,end);
      }
      break;
    case CT_ACO: 
      {
	ACO_Felsenstein fel(order,phylogeny,*Q,A,seqType,numThreads);
	LL=fel.logLikelihood(begin,end);
      }
      break;
    case CT_HOG: 
      {
	HOG_Felsenstein fel(order/2,phylogeny,*Q,A,seqType,nmerTable,
			    numThreads);
	LL=fel.logLikelihood(begin,end);
      }
      break;
    }
  return LL;
}



double NthOrdFelsenstein::logLikelihood3(int begin,int end) 
{
  return logLikelihoodInPhase(begin,end,2);
}



double NthOrdFelsenstein::logLikelihoodInPhase(int begin,int end,int phase) 
{
  double LL=0;
  switch(contextType) 
    {
    case CT_MP: 
      {
	FitchFelsenstein fel(order,phylogeny,*Q,A,seqType,numThreads);
	LL=fel.logLikelihoodInPhase(begin,end,phase);
      }
      break;
    case CT_TRCO: 
      {
	TRCO_Felsenstein fel(order,phylogeny,*Q,A,seqType,numThreads);
	LL=fel.logLikelihoodInPhase(begin,end,phase);
      }
      break;
    case CT_RCO: 
      {
	RCO_Felsenstein fel(order,phylogeny,*Q,A,seqType,numThreads);
	LL=fel.logLikelihoodInPhase(begin,end,phase);
      }
      break;
    case CT_LCO: 
      {
	LCO_Felsenstein fel(order,phylogeny,*Q,A,seqType,numThreads);
	LL=fel.logLikelihoodInPhase(begin,end,phase);
      }
      break;
    case CT_ACO: 
      {
	ACO_Felsenstein fel(order,phylogeny,*Q,A,seqType,numThreads);
	LL=fel.logLikelihoodInPhase(begin,end,phase);
      }
      break;
    case CT_HOG: 
      {
	HOG_Felsenstein fel(order/2,phylogeny,*Q,A,seqType,nmerTable,
			    numThreads);
	LL=fel.logLikelihoodInPhase(begin,end,phase);
      }
      break;
    }
  return LL;
}





