/****************************************************************
 ACO_Felsenstein.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ACO_Felsenstein_H
#define INCL_ACO_Felsenstein_H
#include <iostream>
#include "BOOM/Alphabet.H"
#include "Phylogeny.H"
#include "NthOrdRateMatrix.H"
#include "AlignmentNmerTable.H"
#include "BOOM/MultSeqAlignment.H"
#include "HOG_Felsenstein.H"
using namespace std;
using namespace BOOM;


class ACO_Felsenstein : public HOG_Felsenstein
{
public:
  ACO_Felsenstein(Phylogeny &,NthOrdRateMatrix *,MultSeqAlignment &,
		  AlignmentNmerTable *,int numThreads=1);
  virtual HOG_Likelihooder *getLikelihooder(int numNodes,int numAlpha,
		    Alphabet &alpha,MultSeqAlignment &A,int firstCol,
		    int lastCol,RootNode *root,
		    AlphabetMap &alphabetMap,MolecularSequenceType);
};

#endif

