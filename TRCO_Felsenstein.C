/****************************************************************
 TRCO_Felsenstein.C
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
#include "TRCO_Felsenstein.H"
using namespace std;
using namespace BOOM;


// -----------------------------------------------------------------
// This is where the actual likelihood algorithm is implemented...
// -----------------------------------------------------------------
struct Likelihooder : public TreeVisitor
{
public:
  Array2D<double> L; // the dynamic programming matrix
  MultSeqAlignment &A;
  int numNodes, numAlpha, column, order;
  const Alphabet &alpha;
  double likelihood;
  Sequence rootContext;
  unsigned contextCode, contextLength;
  BitSet gapSymbols;
  AlphabetMap &alphabetMap;
  Array1D<double> V;
  RootNode *root;
  
  Likelihooder(int numNodes,int numAlpha,const Alphabet &alpha,
	       AlphabetMap &aMap,MultSeqAlignment &A,int column,
	       RootNode *root,int order)
    : A(A), numNodes(numNodes), numAlpha(numAlpha),
      alpha(alpha), L(numNodes,numAlpha), column(column),
      gapSymbols(A.getGapSymbols()), alphabetMap(aMap),
      V(numAlpha), order(order), root(root)
  {
    // ctor
    
    reinit(column);
  }
  
  void reinit(int newCol)
  {
    column=newCol;
    int contextBegin=column-order;
    if(contextBegin<0) contextBegin=0;
    int contextLen=column-contextBegin;
    int rootID=root->getID();
    AlignmentSeq &rootTrack=A.getIthTrack(rootID);
    rootContext.clear();
    if(contextLen>0)
      rootTrack.getSeq().getSubsequence(contextBegin,contextLen,
					rootContext);
    int pos=MultSeqAlignment::rightmostGapPos(rootContext,gapSymbols);
    if(pos>=0) {
      Sequence temp;
      int newLen=contextLen-pos-1;
      if(newLen>0)
	rootContext.getSubsequence(pos+1,newLen,temp);
      rootContext=temp;
    }
    contextCode=rootContext.asInt(alpha);
    contextLength=rootContext.getLength();
  }

  virtual void processNode(LeafNode &u) // DO NOT CHANGE THIS FUNCTION
  {
    int id=u.getID();
    //cout<<"leaf: "<<id<<"="<<A.getIthTrack(id).getName()<<endl;
    Symbol unmapped=A.getIthTrack(id)[column];
    Symbol a=alphabetMap(unmapped);
    Array2D<double>::RowIn2DArray<double> row=L[id];
    if(gapSymbols.isMember(unmapped))
      for(Symbol i=0 ; i<numAlpha; ++i)
	row[i]=log(1);// missing data -- same as Seipel & Haussler
    else {
      /*
      cout<<A.getIthTrack(id).isTarget()<<endl;
      cout<<A.getIthTrack(id).getName()<<endl;
      cout<<A.getIthTrack(id).getID()<<endl;
      cout<<int(A.getIthTrack(id)[column])<<" "<<int(unmapped)<<endl;
      throw "bad";
      */
      for(Symbol i=0 ; i<numAlpha; ++i) 
	row[i]=(i==a ? log(1) : log(0));
    }
  }
  
  virtual void processNode(RootNode &u) {
    int id=u.getID();
    //cout<<"root: "<<id<<"="<<A.getIthTrack(id).getName()<<endl;
    Symbol unmapped=A.getIthTrack(id)[column];
    Symbol a=alphabetMap(unmapped);
    int childID=u.getChild()->getID();
    NthOrdSubstMatrix &Pt=*u.getSubstMatrix();
    if(gapSymbols.isMember(unmapped)) {
      for(int a=0 ; a<numAlpha ; ++a)
	V[a]=processInternalChild((Symbol)a,childID,Pt);
      likelihood=sumLogProbs(V);
    }
    else
      likelihood=processInternalChild(a,childID,Pt);
  }
  
  virtual void processNode(InternalNode &u) {
    int id=u.getID();
    //cout<<"internal: "<<id<<"="<<A.getIthTrack(id).getName()<<endl;
    Array2D<double>::RowIn2DArray<double> row=L[id];
    for(Symbol a=0 ; a<numAlpha ; ++a) {
      int left=u.getLeft()->getID(), right=u.getRight()->getID();
      NthOrdSubstMatrix &leftPt=*u.getLeftSubstMatrix();
      NthOrdSubstMatrix &rightPt=*u.getRightSubstMatrix();
      row[a]=
	processInternalChild(a,left,leftPt)+
	processInternalChild(a,right,rightPt);
    }
  }
  
  double processInternalChild(Symbol parentSymbol,int child,
			      NthOrdSubstMatrix &noPt) {
    SubstitutionMatrix &Pt=*noPt.lookup(contextCode,order);
    Array2D<double>::RowIn2DArray<double> row=L[child];
    for(Symbol b=0 ; b<numAlpha ; ++b)
      V[b]=row[b]+Pt(parentSymbol,(Symbol)b);
    return sumLogProbs(V);
  }
  
  double getLikelihood()
  {
    return likelihood;
  }
};



// -----------------------------------------------------------------
// constructor
// -----------------------------------------------------------------
TRCO_Felsenstein::TRCO_Felsenstein(int order,
				     Phylogeny &P,
				     NthOrdRateMatrix &Q,
				     MultSeqAlignment &A,
				     MolecularSequenceType T,
				     int numThreads)
    : phylogeny(P),
      Q(Q),
      alignment(A),
      seqType(T),
      alphabet(Q.getAlphabet()),
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
void TRCO_Felsenstein::attachPhylogenyNodes()
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
double TRCO_Felsenstein::logLikelihood(int column)
{
  PhylogenyNode *root=phylogeny.getRoot();
  if(root->getNodeType()!=ROOT_NODE) {
    cout<<"You must re-root the phylogeny first"<<endl;
    throw 0;
  }
  Likelihooder L(largestID+1,numAlpha,alphabet,Q.getAlphabetMap(),
		 alignment,column,root,order);
  phylogeny.postorderTraversal(L);
  double LL=L.getLikelihood();
  return LL;
}



// -----------------------------------------------------------------
// log-likelihood summed over all columns
// -----------------------------------------------------------------
double TRCO_Felsenstein::logLikelihood()
{
  if(numThreads>1)
    return logLikelihood_threaded();
  
  PhylogenyNode *root=phylogeny.getRoot();
  double LL=0;
  int numColumns=alignment.getLength();
  Likelihooder L(largestID+1,numAlpha,alphabet,Q.getAlphabetMap(),
		 alignment,0,root,order);
  for(int i=0 ; i<numColumns ; ++i) {
    //LL+=logLikelihood(i);
    L.reinit(i);
    phylogeny.postorderTraversal(L);
    LL+=L.getLikelihood();
  }
  return LL;
}



// -----------------------------------------------------------------
// log-likelihood for a range of columns
// -----------------------------------------------------------------
double TRCO_Felsenstein::logLikelihood(int begin,int end)
{
  double L=0;
  for(int i=begin ; i<end ; ++i) {
    //cout<<"column "<<i<<endl;
    L+=logLikelihood(i);
  }
  return L;
}



double TRCO_Felsenstein::logLikelihood3(int begin,int end)
{
  return logLikelihoodInPhase(begin,end,2);
}



double TRCO_Felsenstein::logLikelihoodInPhase(int begin,int end,int phase)
{
  double L=0;
  for(int i=begin ; i<end ; ++i) {
    if(i%3==phase)
      L+=logLikelihood(i);
  }
  return L;
}



double TRCO_Felsenstein::logLikelihood_bestFrame(int begin,int end)
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
// multi-threaded version
// -----------------------------------------------------------------
double TRCO_Felsenstein::logLikelihood_threaded()
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
        FelThread *thread=new FelThread(*this,begin,end);
        threads[i]=thread;
        thread->start();
        nextCol=end;
    }

    // Wait for all the threads to finish
    double L=0;
    for(int i=0 ; i<numThreads ; ++i)
    {
        FelThread *thread=threads[i];
        thread->join();
        L+=thread->getLogLikelihood();
        delete thread;
    }

    // Return the likelihood
    return L;
}



/****************************************************************
                           FelThread methods
 ****************************************************************/
FelThread::FelThread(TRCO_Felsenstein &alg,int begin,int end)
    : alg(alg),
      begin(begin),
      end(end)
{
    // ctor
}



void FelThread::f()
{
    L=0.0;
    for(int i=begin ; i<end ; ++i)
        L+=alg.logLikelihood(i);
}



double FelThread::getLogLikelihood()
{
    return L;
}



