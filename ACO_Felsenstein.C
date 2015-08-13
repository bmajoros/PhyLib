/****************************************************************
 ACO_Felsenstein.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/DnaAlphabet.H"
#include "BOOM/AminoAlphabet.H"
#include "BOOM/Array2D.H"
#include "BOOM/SumLogProbs.H"
#include "ACO_Felsenstein.H"
#include "SafeAdd.H"
#include "DegenerateDnaMatch.H"
#include "BOOM/BitSet.H"
using namespace std;
using namespace BOOM;



// -----------------------------------------------------------------
// This is where the actual likelihood algorithm is implemented...
// -----------------------------------------------------------------
struct ACO_Likelihooder : public TreeVisitor
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
    
  ACO_Likelihooder(int numNodes,int numAlpha,Alphabet &alpha,
		   MultSeqAlignment &A,int firstCol,int lastCol,
		   RootNode *root,AlphabetMap &alphabetMap,
		   MolecularSequenceType seqType)
    : A(A), numNodes(numNodes), numAlpha(numAlpha),alpha(alpha), 
      firstCol(firstCol), lastCol(lastCol), alphabetMap(alphabetMap),
      seqType(seqType), gapSymbols(A.getGapSymbols())
  {
    numCols=lastCol-firstCol+1;
    numNmers=pow((float)alphabetMap.getRangeSize(),(float)numCols);
    L.resize(numNodes,numNmers);
  }
  
  virtual void processNode(LeafNode &u) {
    int id=u.getID();
    Array2D<double>::RowIn2DArray<double> row=L[id];
    Sequence &track=A.getIthTrack(id).getSeq();
    Sequence leafNmer, nmer;
    track.getSubsequence(firstCol,numCols,leafNmer);
    if(MultSeqAlignment::rightmostGapPos(leafNmer,gapSymbols)>=0)
      for(int i=0 ; i<numNmers ; ++i) {
	nmer.fromInt(i,numCols,alphabetMap);
	row[i]=
	  degenerateDnaMatch(nmer,leafNmer,gapSymbols) ? 0 :
	  NEGATIVE_INFINITY;
      }
    else {
      row.setAllTo(NEGATIVE_INFINITY);
      int index=leafNmer.asInt(alphabetMap);
      row[index]=0;
    }
  }
  
  virtual void processNode(RootNode &u) {
    int id=u.getID();
    int childID=u.getChild()->getID();
    NthOrdSubstMatrix &Pt=*u.getSubstMatrix();
    Sequence &track=A.getIthTrack(id).getSeq();
    Sequence rootNmer;
    track.getSubsequence(firstCol,numCols,rootNmer);
    likelihood=processInternalChild(rootNmer,childID,Pt);
    /*
      Vector<double> V;
      for(int i=0 ; i<numNmers ; ++i) {
      nmer.fromInt(i,numCols,alphabetMap);
      if(degenerateDnaMatch(nmer,rootNmer,N))
      V.push_back(processInternalChild(nmer,childID,Pt));
      }            
      likelihood=sumLogProbs(V);//-log(V.size());
    */
  }
  
  virtual void processNode(InternalNode &u) {
    int id=u.getID();
    Array2D<double>::RowIn2DArray<double> row=L[id];
    int left=u.getLeft()->getID(), right=u.getRight()->getID();
    NthOrdSubstMatrix &leftPt=*u.getLeftSubstMatrix();
    NthOrdSubstMatrix &rightPt=*u.getRightSubstMatrix();
    Sequence nmer;
    for(int i=0 ; i<numNmers ; ++i) {
      nmer.fromInt(i,numCols,alphabetMap);
      row[i]=
	processInternalChild(nmer,left,leftPt)+
	processInternalChild(nmer,right,rightPt);
    }
  }
  
  double processInternalChild(Sequence &parentNmer,int childID,
			      NthOrdSubstMatrix &Pt) {
    Array2D<double>::RowIn2DArray<double> row=L[childID];
    Array1D<double> V(numNmers);
    Sequence childNmer;
    for(int i=0 ; i<numNmers ; ++i) {
      childNmer.fromInt(i,numCols,alphabetMap);
      double joint=jointProb(parentNmer,childNmer,Pt);
      V[i]=safeAdd(joint,row[i]);
    }
    double ll=sumLogProbs(V);
    return ll;
  }
  
  // jointProb() precodition: neither nmer is degenerate!
  double jointProb(Sequence &parentNmer,Sequence &childNmer,
		   NthOrdSubstMatrix &noPt) {
    double logP=0;
    Sequence parentContext, childContext;
    for(int i=0 ; i<numCols ; ++i) {
      Symbol parentSymbol=parentNmer[i], childSymbol=childNmer[i];
      parentNmer.getSubsequence(0,i,parentContext);
      childNmer.getSubsequence(0,i,childContext);
      SubstitutionMatrix &Pt=*noPt.lookup(parentContext);
      logP+=Pt(parentSymbol,childSymbol);
    }
    return logP;
  }
  
  double getLikelihood() {return likelihood;}
};



// -----------------------------------------------------------------
// constructor
// -----------------------------------------------------------------
ACO_Felsenstein::ACO_Felsenstein(int order,
				     Phylogeny &P,
				     NthOrdRateMatrix &Q,
				     MultSeqAlignment &A,
				     MolecularSequenceType T,
				     int numThreads)
    : phylogeny(P),
      Q(Q),
      alignment(A),
      seqType(T),
      alphabet(T==DNA ? (Alphabet&) DnaAlphabet::global() : 
               (Alphabet&) AminoAlphabet::global()),
      numThreads(numThreads),
      threads(numThreads),
      order(order),
      alphabetMap(Q.getAlphabetMap())
{
    numAlpha=alphabet.size();
    attachPhylogenyNodes();
}



// -----------------------------------------------------------------
// attach phylogeny nodes
// -----------------------------------------------------------------
void ACO_Felsenstein::attachPhylogenyNodes()
{
    struct Attacher : public TreeVisitor
    {
        MultSeqAlignment &A;
        NthOrdRateMatrix &Q;
        int numNodes, largestID;
        
        Attacher(MultSeqAlignment &A,NthOrdRateMatrix &Q) 
            : A(A), Q(Q), numNodes(0), largestID(0) {}
        virtual void processNode(InternalNode &V) 
            {
                V.setLeftSubstMatrix(Q.instantiate(V.getLeftDistance()));
                V.setRightSubstMatrix(Q.instantiate(V.getRightDistance()));
                int id=V.getID();
                if(id>largestID) largestID=id;
                ++numNodes;
            }
        virtual void processNode(LeafNode &V) 
            {
                int id=V.getID();
                if(id>largestID) largestID=id;
                ++numNodes;
            }
        virtual void processNode(RootNode &V)
            { 
                V.setSubstMatrix(Q.instantiate(V.getBranchLength()));
                int id=V.getID();
                if(id>largestID) largestID=id;
                ++numNodes;
            }
    };
    
    Attacher A(alignment,Q);
    phylogeny.postorderTraversal(A);
    numNodes=A.numNodes;
    largestID=A.largestID;
}



// -----------------------------------------------------------------
// (conditional) log-likelihood for a single column
// -----------------------------------------------------------------
double ACO_Felsenstein::logLikelihood(int column)
{
  // Make sure the phylogeny has been properly re-rooted
  PhylogenyNode *root=phylogeny.getRoot();
  if(root->getNodeType()!=ROOT_NODE) {
    cout<<"You must re-root the phylogeny first"<<endl;
    throw 0;
  }
  int rootID=root->getID();

  // Compute dimensions of the slices for joint probability computation
  int sliceBegin=column-order;
  if(sliceBegin<0) sliceBegin=0;
  int bigSliceEnd=column, smallSliceEnd=column-1;
  int bigSliceWidth=bigSliceEnd-sliceBegin+1;
  int smallSliceWidth=smallSliceEnd-sliceBegin+1;
  
  // Compute joint probability for N columns
  ACO_Likelihooder bigSlice(largestID+1,numAlpha,alphabet,alignment,
                            sliceBegin,bigSliceEnd,root,alphabetMap,seqType);
  phylogeny.postorderTraversal(bigSlice);
  double bigJoint=bigSlice.getLikelihood();

  if(smallSliceWidth<1) return bigJoint; // equals 0th-order conditional
  else {
      // Compute joint probability for N-1 columns
      ACO_Likelihooder smallSlice(largestID+1,numAlpha,alphabet,alignment,
                                  sliceBegin,smallSliceEnd,root,alphabetMap,
				  seqType);
      phylogeny.postorderTraversal(smallSlice);
      double smallJoint=smallSlice.getLikelihood();
      
      // Form the conditional probability by dividing the joint for the big
      // slice by the joint for the smaller slice:
      return bigJoint-smallJoint;
  }
}



// -----------------------------------------------------------------
// log-likelihood summed over all columns
// -----------------------------------------------------------------
double ACO_Felsenstein::logLikelihood()
{
  //return logLikelihood_threaded();
  
  double L=0;
  int numColumns=alignment.getLength();
  for(int i=0 ; i<numColumns ; ++i)
    L+=logLikelihood(i);
  return L;
}


double ACO_Felsenstein::logLikelihood3(int begin,int end)
{
  return logLikelihoodInPhase(begin,end,2);
}


double ACO_Felsenstein::logLikelihoodInPhase(int begin,int end,int phase)
{
  double L=0;
  for(int i=begin ; i<end ; ++i)
    if(i%3==phase)
      L+=logLikelihood(i);
  return L;
}


double ACO_Felsenstein::logLikelihood_bestFrame(int begin,int end)
{
  double bestLL=NEGATIVE_INFINITY;
  int bestSampleSize=0;
  for(int frame=0 ; frame<3 ; ++frame) {
    double frameLL=0;
    int frameSampleSize=0;
    for(int i=begin ; i<end ; ++i)
      if(i%3==frame) {
	frameLL+=logLikelihood(i);
	++frameSampleSize;
      }
    if(bestSampleSize==0 || 
       frameLL/frameSampleSize>bestLL/bestSampleSize) {
      bestLL=frameLL;
      bestSampleSize=frameSampleSize;
    }
  }
  return bestLL;
}



// -----------------------------------------------------------------
// log-likelihood for a range of columns
// -----------------------------------------------------------------
double ACO_Felsenstein::logLikelihood(int begin,int end)
{
  double L=0;
  for(int i=begin ; i<end ; ++i)
    L+=logLikelihood(i);
  return L;
}



// -----------------------------------------------------------------
// multi-threaded version
// -----------------------------------------------------------------
double ACO_Felsenstein::logLikelihood_threaded()
{
    // Spawn all the threads
    int numColumns=alignment.getLength();
    int colsPerThread=numColumns/numThreads;
    int nextCol=0;
    for(int i=0 ; i<numThreads ; ++i)
    {
        int begin=nextCol;
        int end=begin+colsPerThread;
        if(i+1==numThreads) end=numColumns;
        ACO_FelThread *thread=new ACO_FelThread(*this,begin,end);
        threads[i]=thread;
        thread->start();
        nextCol=end;
    }

    // Wait for all the threads to finish
    double L=0;
    for(int i=0 ; i<numThreads ; ++i)
    {
        ACO_FelThread *thread=threads[i];
        thread->join();
        L+=thread->getLogLikelihood();
        delete thread;
    }

    // Return the likelihood
    return L;
}



/****************************************************************
                           ACO_FelThread methods
 ****************************************************************/
ACO_FelThread::ACO_FelThread(ACO_Felsenstein &alg,int begin,int end)
    : alg(alg),
      begin(begin),
      end(end)
{
    // ctor
}



void ACO_FelThread::f()
{
    L=0.0;
    for(int i=begin ; i<end ; ++i)
        L+=alg.logLikelihood(i);
}



double ACO_FelThread::getLogLikelihood()
{
    return L;
}



