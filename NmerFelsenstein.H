/****************************************************************
 NmerFelsenstein.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_NmerFelsenstein_H
#define INCL_NmerFelsenstein_H
#include <iostream>
#include "BOOM/Alphabet.H"
#include "Phylogeny.H"
#include "NmerRateMatrix.H"
#include "AlignmentNmerTable.H"
#include "BOOM/MultSeqAlignment.H"
#include "HOG_Felsenstein.H"
using namespace std;
using namespace BOOM;


class NmerFelsenstein : public HOG_Felsenstein
{
public:
  NmerFelsenstein(Phylogeny &,NmerRateMatrix *,MultSeqAlignment &,
		  AlignmentNmerTable *,int numThreads=1);
  virtual HOG_Likelihooder *getLikelihooder(int numNodes,int numAlpha,
		    Alphabet &alpha,MultSeqAlignment &A,int firstCol,
		    int lastCol,RootNode *root,
		    AlphabetMap &alphabetMap,MolecularSequenceType);
};

#endif

