/****************************************************************
 HOG_Felsenstein.C
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
#include "BOOM/BitSet.H"
#include "HOG_Felsenstein.H"
#include "SafeAdd.H"
#include "DegenerateDnaMatch.H"
using namespace std;
using namespace BOOM;



// -----------------------------------------------------------------
// HOG_Likelihooder methods
// -----------------------------------------------------------------
    
HOG_Likelihooder::HOG_Likelihooder(int numNodes,int numAlpha,Alphabet &alpha,
				   MultSeqAlignment &A,int firstCol,
				   int lastCol,RootNode *root,
				   AlphabetMap &alphabetMap,
				   MolecularSequenceType seqType,
				   AlignmentNmerTable *nmerTable)
  : A(A), numNodes(numNodes), numAlpha(numAlpha), alpha(alpha), 
    firstCol(firstCol), lastCol(lastCol), alphabetMap(alphabetMap), 
    seqType(seqType), gapSymbols(A.getGapSymbols()),
    nmerTable(nmerTable)
{
  numCols=lastCol-firstCol+1;
  numNmers=pow((float)alphabetMap.getRangeSize(),(float)numCols);
  L.resize(numNodes,numNmers);
}



void HOG_Likelihooder::reinit(int first,int last)
{
  firstCol=first;
  lastCol=last;
  numCols=lastCol-firstCol+1;
}
    


void HOG_Likelihooder::processNode(LeafNode &u) 
{
  int id=u.getID();
  Array2D<double>::RowIn2DArray<double> row=L[id];
  Sequence &track=A.getIthTrack(id).getSeq();
  Sequence leafNmer, nmer;//###could allocate in reinit()
  track.getSubsequence(firstCol,numCols,leafNmer);
  for(int i=0 ; i<numNmers ; ++i) {
    nmer.fromInt(i,numCols,alphabetMap);
    row[i]=
      degenerateDnaMatch(nmer,leafNmer,gapSymbols) ? 0 :
      NEGATIVE_INFINITY;
  }
}


void HOG_Likelihooder::processNode(RootNode &u) 
{
  int id=u.getID();
  int childID=u.getChild()->getID();
  NthOrdSubstMatrix &Pt=*u.getSubstMatrix();
  Sequence &track=A.getIthTrack(id).getSeq();
  Sequence rootNmer;
  track.getSubsequence(firstCol,numCols,rootNmer);
  likelihood=processInternalChild(rootNmer,childID,Pt);
}



void HOG_Likelihooder::processNode(InternalNode &u) 
{
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



double HOG_Likelihooder::processInternalChild(Sequence &parentNmer,
					      int childID,
					      NthOrdSubstMatrix &Pt) 
{
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



// jointProb() precondition: neither nmer is degenerate!
double HOG_Likelihooder::jointProb(Sequence &parentNmer,Sequence &childNmer,
				   NthOrdSubstMatrix &dualPt) 
{
  double logP=0;
  Sequence parentContext, childContext;
  for(int i=0 ; i<numCols ; ++i) {
    Symbol parentSymbol=parentNmer[i], childSymbol=childNmer[i];
    parentNmer.getSubsequence(0,i,parentContext);
    childNmer.getSubsequence(0,i,childContext);
    SubstitutionMatrix &Pt=*dualPt.lookup(parentContext,childContext);
    logP+=log(Pt(parentSymbol,childSymbol));
  }
  return logP;
}



double HOG_Likelihooder::getLikelihood() 
{
  return likelihood;
}



// -----------------------------------------------------------------
// HOG_Felsenstein methods
// -----------------------------------------------------------------

HOG_Likelihooder *HOG_Felsenstein::getLikelihooder(int numNodes,int numAlpha,
		    Alphabet &alpha,MultSeqAlignment &A,int firstCol,
		    int lastCol,RootNode *root,
		   AlphabetMap &alphabetMap,MolecularSequenceType seqType)
{
  return new HOG_Likelihooder(numNodes,numAlpha,alpha,A,firstCol,lastCol,
			      root,alphabetMap,seqType,nmerTable);
}



HOG_Felsenstein::HOG_Felsenstein(int order,
				 Phylogeny &P,
				 NthOrdRateMatrix &Q,
				 MultSeqAlignment &A,
				 MolecularSequenceType T,
				 AlignmentNmerTable *nmerTable,
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
      alphabetMap(Q.getAlphabetMap()),
      gapSymbols(A.getGapSymbols()),
      nmerTable(nmerTable)
{
  numAlpha=alphabet.size();

  // Make sure the phylogeny has been properly re-rooted
  root=phylogeny.getRoot();
  if(root->getNodeType()!=ROOT_NODE) {
    cout<<"You must re-root the phylogeny first"<<endl;
    throw 0;
  }
  
  int rootID=root->getID();
  attachPhylogenyNodes();
  rootTrack=&alignment.getIthTrack(rootID);
}



// -----------------------------------------------------------------
// attach phylogeny nodes
// -----------------------------------------------------------------
void HOG_Felsenstein::attachPhylogenyNodes()
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
double HOG_Felsenstein::logLikelihood(int column)
{
  // Compute dimensions of the slices for joint probability computation
  int sliceBegin=column-order;
  if(sliceBegin<0) sliceBegin=0;
  int bigSliceEnd=column, smallSliceEnd=column-1;
  int bigSliceWidth=bigSliceEnd-sliceBegin+1;

  // Make sure there are no gaps in the root sequence for this slice
  /*
  Sequence rootNmer;
  rootTrack->getSeq().getSubsequence(sliceBegin,bigSliceWidth,rootNmer);
  int pos=MultSeqAlignment::rightmostGapPos(rootNmer,gapSymbols);
  if(pos==bigSliceWidth-1) return 0;
  if(pos>=0) {
    sliceBegin+=(pos+1);
    bigSliceWidth=bigSliceEnd-sliceBegin+1;
  }
  */
  
  // Compute joint probability for N columns
  bigSlice->reinit(sliceBegin,bigSliceEnd);
  phylogeny.postorderTraversal(*bigSlice);
  double bigJoint=bigSlice->getLikelihood();

  int smallSliceWidth=smallSliceEnd-sliceBegin+1;
  if(smallSliceWidth<1) return bigJoint; // equals 0th-order conditional
  else {
      // Compute joint probability for N-1 columns
      smallSlice->reinit(sliceBegin,smallSliceEnd);
      phylogeny.postorderTraversal(*smallSlice);
      double smallJoint=smallSlice->getLikelihood();

      // Form the conditional probability by dividing the joint for the big
      // slice by the joint for the smaller slice:
      return bigJoint-smallJoint;
  }
}



// -----------------------------------------------------------------
// log-likelihood summed over all columns
// -----------------------------------------------------------------
double HOG_Felsenstein::logLikelihood()
{
  bigSlice=getLikelihooder(largestID+1,numAlpha,alphabet,alignment,0,order,
			   root,alphabetMap,seqType);
  smallSlice=getLikelihooder(largestID+1,numAlpha,alphabet,alignment,0,
			     order-1, // looks wrong, but it's really correct
			     root,alphabetMap,seqType);

  if(numThreads>1)
    return logLikelihood_threaded();
  
  double L=0;
  int numColumns=alignment.getLength();
  for(int i=0 ; i<numColumns ; ++i)
    L+=logLikelihood(i);

  delete bigSlice;
  delete smallSlice;
  return L;
}



// -----------------------------------------------------------------
// log-likelihood for a range of columns
// -----------------------------------------------------------------
double HOG_Felsenstein::logLikelihood(int begin,int end)
{
  bigSlice=getLikelihooder(largestID+1,numAlpha,alphabet,alignment,0,order,
			   root,alphabetMap,seqType);
  smallSlice=getLikelihooder(largestID+1,numAlpha,alphabet,alignment,0,
			     order-1, // looks wrong, but it's actually right
			     root,alphabetMap,seqType);
  double L=0;
  for(int i=begin ; i<end ; ++i)
    L+=logLikelihood(i);
  delete bigSlice;
  delete smallSlice;
  return L;
}



double HOG_Felsenstein::logLikelihood3(int begin,int end)
{
  return logLikelihoodInPhase(begin,end,2);
}



double HOG_Felsenstein::logLikelihoodInPhase(int begin,int end,int phase)
{
  bigSlice=getLikelihooder(largestID+1,numAlpha,alphabet,alignment,0,order,
			   root,alphabetMap,seqType);
  smallSlice=getLikelihooder(largestID+1,numAlpha,alphabet,alignment,0,order-1,
			     root,alphabetMap,seqType);
  double L=0;
  for(int i=begin ; i<end ; ++i)
    if(i%3==phase)
      L+=logLikelihood(i);
  delete bigSlice;
  delete smallSlice;
  return L;
}



double HOG_Felsenstein::logLikelihood_bestFrame(int begin,int end)
{
  bigSlice=getLikelihooder(largestID+1,numAlpha,alphabet,alignment,0,order,
			   root,alphabetMap,seqType);
  smallSlice=getLikelihooder(largestID+1,numAlpha,alphabet,alignment,0,order-1,
			     root,alphabetMap,seqType);
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
  delete bigSlice;
  delete smallSlice;
  return bestLL;
}


// -----------------------------------------------------------------
// multi-threaded version
// -----------------------------------------------------------------
double HOG_Felsenstein::logLikelihood_threaded()
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
        HOG_FelThread *thread=new HOG_FelThread(*this,begin,end);
        threads[i]=thread;
        thread->start();
        nextCol=end;
    }

    // Wait for all the threads to finish
    double L=0;
    for(int i=0 ; i<numThreads ; ++i)
    {
        HOG_FelThread *thread=threads[i];
        thread->join();
        L+=thread->getLogLikelihood();
        delete thread;
    }

    // Return the likelihood
    return L;
}



/****************************************************************
                           HOG_FelThread methods
 ****************************************************************/
HOG_FelThread::HOG_FelThread(HOG_Felsenstein &alg,int begin,int end)
    : alg(alg),
      begin(begin),
      end(end)
{
    // ctor
}



void HOG_FelThread::f()
{
    L=0.0;
    for(int i=begin ; i<end ; ++i)
        L+=alg.logLikelihood(i);
}



double HOG_FelThread::getLogLikelihood()
{
    return L;
}



