/****************************************************************
 FitchFelsenstein.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/DnaAlphabet.H"
#include "BOOM/AminoAlphabet.H"
#include "BOOM/Array2D.H"
#include "BOOM/Vector.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/BitSet.H"
#include "FitchFelsenstein.H"
#include "SafeAdd.H"
using namespace std;
using namespace BOOM;


// -----------------------------------------------------------------
// This is where the actual likelihood algorithm is implemented...
// -----------------------------------------------------------------
struct FitchLikelihooder : public TreeVisitor
{
public:
  Array2D<double> L; // the dynamic programming matrix
  MultSeqAlignment &A;
  int numNodes, numAlpha, column, order;
  Alphabet &alpha;
  double likelihood;
  RootNode *root; // ### for debugging only
  BitSet gapSymbols;
  Array1D<double> V;
    
  FitchLikelihooder(int numNodes,int numAlpha,Alphabet &alpha,
		    MultSeqAlignment &A,int column,RootNode *root,int order)
    : A(A), numNodes(numNodes), numAlpha(numAlpha),
      alpha(alpha), L(numNodes,numAlpha), column(column),
      order(order), root(root), likelihood(NEGATIVE_INFINITY),
      gapSymbols(A.getGapSymbols()), V(numAlpha)
  {
    L.setAllTo(NEGATIVE_INFINITY);
    V.setAllTo(NEGATIVE_INFINITY);
  }
  
  void reinit(int col)
  {
    L.setAllTo(NEGATIVE_INFINITY);
    column=col;
  }

  virtual void processNode(LeafNode &u) // DO NOT CHANGE THIS FUNCTION
  {
    int id=u.getID();
    Symbol a=A.getIthTrack(id)[column];
    Array2D<double>::RowIn2DArray<double> row=L[id];
    if(gapSymbols.isMember(a))
      for(Symbol i=0 ; i<numAlpha; ++i)
	row[i]=0;// missing data -- same as Seipel & Haussler
    else
      for(Symbol i=0 ; i<numAlpha; ++i) 
	row[i]=(i==a ? 0 : NEGATIVE_INFINITY);
  }
  
  virtual void processNode(RootNode &u) {
    // PRECONDITION: no gaps in target sequence!
    int id=u.getID();
    Symbol a=A.getIthTrack(id)[column];
    int childID=u.getChild()->getID();
    NthOrdSubstMatrix &Pt=*u.getSubstMatrix();
    likelihood=processInternalChild(a,id,childID,Pt);
  }
  
  virtual void processNode(InternalNode &u) {
    int id=u.getID();
    NthOrdSubstMatrix &leftPt=*u.getLeftSubstMatrix();
    NthOrdSubstMatrix &rightPt=*u.getRightSubstMatrix();
    int left=u.getLeft()->getID(), right=u.getRight()->getID();
    Array2D<double>::RowIn2DArray<double> row=L[id];
    Symbol a=A.getIthTrack(id)[column];
    if(gapSymbols.isMember(a))
      for(a=0 ; a<numAlpha ; ++a) {
	if(gapSymbols.isMember(a)) continue;
	row[a]=
	  processInternalChild(a,id,left,leftPt)+
	  processInternalChild(a,id,right,rightPt);
      }
    else {
      for(Symbol b=0 ; b<numAlpha ; ++b) row[b]=NEGATIVE_INFINITY;
      double l=processInternalChild(a,id,left,leftPt);
      double r=processInternalChild(a,id,right,rightPt);
      row[a]=l+r;
    }
  }
  
  double processInternalChild(Symbol parentSymbol,int parentID,int child,
			      NthOrdSubstMatrix &noPt) {
    AlignmentSeq &track=A.getIthTrack(parentID);
    AlignmentSeq &childTrack=A.getIthTrack(child);
    int contextBegin=column-order;
    if(contextBegin<0) contextBegin=0;
    int contextLen=column-contextBegin;
    Sequence context;
    track.getSeq().getSubsequence(contextBegin,contextLen,context);//ACO
    //A.getIthTrack(root->getID()).getSeq().getSubsequence(contextBegin,contextLen,context);//TRCO
    int pos=MultSeqAlignment::rightmostGapPos(context,gapSymbols);
    if(pos>=0) {
      throw "this should not happen (FitchFelsenstein)::processInternalChild";
      Sequence temp;
      context.getSubsequence(pos+1,contextLen-pos-1,temp);
      context=temp;
    }
    SubstitutionMatrix &Pt=*noPt.lookup(context,0,-1);
    Array2D<double>::RowIn2DArray<double> row=L[child];
    Symbol childSym=childTrack.getSeq()[column];
    double ll;
    if(gapSymbols.isMember(childSym)){
      V.setAllTo(NEGATIVE_INFINITY);
      for(Symbol a=0 ; a<numAlpha ; ++a) {
	if(!gapSymbols.isMember(a)) {
	  V[(int)a]=safeAdd(row[a],Pt(parentSymbol,a));
	}
      }
      ll=sumLogProbs(V);
    }
    else {
      ll=safeAdd(row[childSym],Pt(parentSymbol,childSym));
    }
    return ll;
  }
  
  double getLikelihood()
  {
    return likelihood;
  }
};



// -----------------------------------------------------------------
// constructor
// -----------------------------------------------------------------
FitchFelsenstein::FitchFelsenstein(int order,
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
      order(order)
{
    numAlpha=alphabet.size();
    attachPhylogenyNodes();
}



// -----------------------------------------------------------------
// attach phylogeny nodes
// -----------------------------------------------------------------
void FitchFelsenstein::attachPhylogenyNodes()
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
// log-likelihood for a single column
// -----------------------------------------------------------------
double FitchFelsenstein::logLikelihood(int column)
{
  PhylogenyNode *root=phylogeny.getRoot();
  if(root->getNodeType()!=ROOT_NODE) {
    cout<<"You must re-root the phylogeny first"<<endl;
    throw 0;
  }
  //int rootID=root->getID();
  FitchLikelihooder L(largestID+1,numAlpha,alphabet,alignment,column,
		      root,order);
  phylogeny.postorderTraversal(L);
  return L.getLikelihood();
}



// -----------------------------------------------------------------
// log-likelihood summed over all columns
// -----------------------------------------------------------------
double FitchFelsenstein::logLikelihood()
{
  FitchLikelihooder L(largestID+1,numAlpha,alphabet,alignment,0,
		      phylogeny.getRoot(),order);
  double LL=0;
  int numColumns=alignment.getLength();
  for(int i=0 ; i<numColumns ; ++i) {
    L.reinit(i);
    phylogeny.postorderTraversal(L);
    LL+=L.getLikelihood();
  }
  return LL;
}



// -----------------------------------------------------------------
// log-likelihood for a range of columns
// -----------------------------------------------------------------
double FitchFelsenstein::logLikelihood(int begin,int end)
{
  FitchLikelihooder L(largestID+1,numAlpha,alphabet,alignment,begin,
		      phylogeny.getRoot(),order);
  double LL=0;
  for(int i=begin ; i<end ; ++i) {
    L.reinit(i);
    phylogeny.postorderTraversal(L);
    LL+=L.getLikelihood();
  }
  return LL;
}



double FitchFelsenstein::logLikelihood3(int begin,int end)
{
  return logLikelihoodInPhase(begin,end,2);
}



double FitchFelsenstein::logLikelihoodInPhase(int begin,int end,int phase)
{
  FitchLikelihooder L(largestID+1,numAlpha,alphabet,alignment,begin,
		      phylogeny.getRoot(),order);
  double LL=0;
  for(int i=begin ; i<end ; ++i)
    if(i%3==phase) {
      L.reinit(i);
      phylogeny.postorderTraversal(L);
      LL+=L.getLikelihood();
    }
  return LL;
}



double FitchFelsenstein::logLikelihood_bestFrame(int begin,int end)
{
  FitchLikelihooder L(largestID+1,numAlpha,alphabet,alignment,begin,
		      phylogeny.getRoot(),order);
  double bestLL=NEGATIVE_INFINITY;
  int bestSampleSize=0;
  for(int frame=0 ; frame<3 ; ++frame) {
    double frameLL=0;
    int frameSampleSize=0;
    for(int i=begin ; i<end ; ++i)
      if(i%3==frame) {
	L.reinit(i);
	phylogeny.postorderTraversal(L);
	frameLL+=L.getLikelihood();
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
// multi-threaded version
// -----------------------------------------------------------------
double FitchFelsenstein::logLikelihood_threaded()
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
        FitchThread *thread=new FitchThread(*this,begin,end);
        threads[i]=thread;
        thread->start();
        nextCol=end;
    }

    // Wait for all the threads to finish
    double L=0;
    for(int i=0 ; i<numThreads ; ++i)
    {
        FitchThread *thread=threads[i];
        thread->join();
        L+=thread->getLogLikelihood();
        delete thread;
    }

    // Return the likelihood
    return L;
}



/****************************************************************
                           FitchThread methods
 ****************************************************************/
FitchThread::FitchThread(FitchFelsenstein &alg,int begin,int end)
    : alg(alg),
      begin(begin),
      end(end)
{
    // ctor
}



void FitchThread::f()
{
    L=0.0;
    for(int i=begin ; i<end ; ++i)
        L+=alg.logLikelihood(i);
}



double FitchThread::getLogLikelihood()
{
    return L;
}



