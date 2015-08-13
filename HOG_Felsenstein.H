/****************************************************************
 HOG_Felsenstein.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_HOG_Felsenstein_H
#define INCL_HOG_Felsenstein_H
#include <iostream>
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Alphabet.H"
#include "BOOM/MolecularSequenceType.H"
#include "BOOM/Vector.H"
#include "BOOM/Thread.H"
#include "BOOM/Array2D.H"
#include "BOOM/BitSet.H"
#include "Phylogeny.H"
#include "NthOrdRateMatrix.H"
#include "FelsensteinInterface.H"
#include "DegenerateDnaMatch.H"
#include "AlignmentNmerTable.H"
using namespace std;
using namespace BOOM;


class HOG_Felsenstein;
class HOG_Likelihooder;


class HOG_FelThread : public Thread
{
    double L;
    int begin, end; // columns
    HOG_Felsenstein &alg;
public:
    HOG_FelThread(HOG_Felsenstein &,int begin,int end);
    virtual void f();
    double getLogLikelihood();
};



class HOG_Felsenstein : public FelsensteinInterface
{
public:
  HOG_Felsenstein(int order,Phylogeny &,NthOrdRateMatrix &,MultSeqAlignment&,
		  MolecularSequenceType,AlignmentNmerTable *,int numThreads=1);
  virtual double logLikelihood(int column);
  virtual double logLikelihood();
  virtual double logLikelihood(int begin,int end);
  virtual double logLikelihood_bestFrame(int begin,int end);
  virtual double logLikelihood3(int begin,int end);//every third position only
  virtual double logLikelihoodInPhase(int begin,int end,int phase);
  virtual HOG_Likelihooder *getLikelihooder(int numNodes,int numAlpha,
		    Alphabet &alpha,MultSeqAlignment &A,int firstCol,
		    int lastCol,RootNode *root,
		    AlphabetMap &alphabetMap,MolecularSequenceType);
protected:
  Phylogeny &phylogeny;
  NthOrdRateMatrix &Q;
  MultSeqAlignment &alignment;
  MolecularSequenceType seqType;
  Alphabet &alphabet;
  int numNodes;
  int numAlpha;
  int largestID;
  int numThreads;
  Array1D<HOG_FelThread*> threads;
  int order;
  AlphabetMap &alphabetMap;
  AlignmentSeq *rootTrack;
  PhylogenyNode *root;
  BitSet gapSymbols;
  HOG_Likelihooder *bigSlice, *smallSlice;
  AlignmentNmerTable *nmerTable;

  void attachPhylogenyNodes();
  double logLikelihood_threaded();
};



// -----------------------------------------------------------------
// This is where the actual likelihood algorithm is implemented...
// -----------------------------------------------------------------
struct HOG_Likelihooder : public TreeVisitor
{
public:
  Array2D<double> L; // the dynamic programming matrix
  MultSeqAlignment &A;
  int numNodes, numAlpha, firstCol, lastCol, numCols, numNmers;
  Alphabet &alpha;
  double likelihood;
  AlphabetMap &alphabetMap;
  MolecularSequenceType seqType;
  BitSet gapSymbols;
  AlignmentNmerTable *nmerTable;

  HOG_Likelihooder(int numNodes,int numAlpha,Alphabet &alpha,
		   MultSeqAlignment &A,int firstCol,int lastCol,
		   RootNode *root,AlphabetMap &alphabetMap,
		   MolecularSequenceType seqType,
		   AlignmentNmerTable *);
  virtual void processNode(LeafNode &u);
  virtual void processNode(RootNode &u);
  virtual void processNode(InternalNode &u);
  double processInternalChild(Sequence &parentNmer,int childID,
			      NthOrdSubstMatrix &Pt);
  virtual double jointProb(Sequence &parentNmer,Sequence &childNmer,
			   NthOrdSubstMatrix &dualPt);
  double getLikelihood();
  void reinit(int firstCol,int lastCol);
};


#endif

