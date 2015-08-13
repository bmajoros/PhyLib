/****************************************************************
 LCO_Felsenstein.C
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
#include "BOOM/GapPatternAlphabet.H"
#include "LCO_Felsenstein.H"
using namespace std;
using namespace BOOM;


// -----------------------------------------------------------------
// This is where the actual likelihood algorithm is implemented...
// -----------------------------------------------------------------
struct LCO_Likelihooder : public TreeVisitor
{
public:
  Array2D<double> L; // the dynamic programming matrix
  MultSeqAlignment &A;
  int numNodes, numAlpha, column;
  const Alphabet &alpha;
  AlphabetMap &alphabetMap;
  double likelihood;
  Sequence noContext;
  BitSet gapSymbols;
  int order;
  RootNode *root;

  LCO_Likelihooder(int numNodes,int numAlpha,const Alphabet &alpha,
		   MultSeqAlignment &A,int column,RootNode *root,
		   int order,AlphabetMap &alphabetMap,
		   const BitSet &gapSymbols)
    : A(A), numNodes(numNodes), numAlpha(numAlpha),
      alpha(alpha), L(numNodes,numAlpha), column(column),
      alphabetMap(alphabetMap), gapSymbols(gapSymbols), order(order),
      root(root)
  {
    // ctor
  }
  
  virtual void processNode(LeafNode &u) // DO NOT CHANGE THIS FUNCTION
  {
    int id=u.getID();
    Symbol unmapped=A.getIthTrack(id)[column];
    Symbol a=alphabetMap(unmapped);
    Array2D<double>::RowIn2DArray<double> row=L[id];
    if(gapSymbols.isMember(unmapped))
      for(Symbol i=0 ; i<numAlpha; ++i)
	row[i]=0; //=log(1); missing data -- same as Seipel & Haussler
    else
      for(Symbol i=0 ; i<numAlpha; ++i) 
	row[i]=(i==a ? 0 : NEGATIVE_INFINITY);
  }
  
  virtual void processNode(RootNode &u) {
    int id=u.getID();
    Symbol a=alphabetMap(A.getIthTrack(id)[column]);
    NthOrdSubstMatrix &Pt=*u.getSubstMatrix();
    likelihood=processInternalChild(a,u.getChild(),Pt,&u);
  }
  
  virtual void processNode(InternalNode &u) {
    int id=u.getID();
    PhylogenyNode *left=u.getLeft(), *right=u.getRight();
    NthOrdSubstMatrix &leftPt=*u.getLeftSubstMatrix();
    NthOrdSubstMatrix &rightPt=*u.getRightSubstMatrix();
    Array2D<double>::RowIn2DArray<double> row=L[id];
    for(Symbol a=0 ; a<numAlpha ; ++a) {
      row[a]=
	processInternalChild(a,left,leftPt,&u)+
	processInternalChild(a,right,rightPt,&u);
    }
  }
  
  double processInternalChild(Symbol parentSymbol,PhylogenyNode *child,
			      NthOrdSubstMatrix &noPt,PhylogenyNode *parent) 
  {
    int contextCode=0, contextLength=0;
    if(child->getNodeType()==LEAF_NODE)
      getContext(child,contextCode,contextLength);
    else //if(parent->getNodeType()==ROOT_NODE)
      getContext(root,contextCode,contextLength);
    SubstitutionMatrix &Pt=*noPt.lookup(contextCode,contextLength);
    Array2D<double>::RowIn2DArray<double> row=L[child->getID()];
    Vector<double> V;
    for(Symbol b=0 ; b<numAlpha ; ++b)
      V.push_back(row[b]+Pt(parentSymbol,b)); // Pt already in log space
    double ll=sumLogProbs(V);
    return ll;
  }
  
  void getContext(PhylogenyNode *node,int &contextCode,int &contextLength)
  {
    int ID=node->getID();
    AlignmentSeq &track=A.getIthTrack(ID);
    int contextBegin=column-order;
    if(contextBegin<0) contextBegin=0;
    contextLength=column-contextBegin;
    if(contextLength==0) {contextCode=0;return;}
    Sequence context;
    track.getSeq().getSubsequence(contextBegin,contextLength,context);
    int pos=MultSeqAlignment::rightmostGapPos(context,gapSymbols);
    if(pos>=0) {
      Sequence temp;
      int begin=pos+1, len=contextLength-pos-1;
      if(begin<order && len>0)
	context.getSubsequence(pos+1,contextLength-pos-1,temp);
      context=temp;
      contextLength=context.getLength();
    }
    contextCode=context.asInt(*alphabetMap.getRange());
  }

  double getLikelihood()
  {
    return likelihood;
  }
};



// -----------------------------------------------------------------
// constructor
// -----------------------------------------------------------------
LCO_Felsenstein::LCO_Felsenstein(int order,
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
void LCO_Felsenstein::attachPhylogenyNodes()
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
double LCO_Felsenstein::logLikelihood(int column)
{
  PhylogenyNode *root=phylogeny.getRoot();
  if(root->getNodeType()!=ROOT_NODE) {
    cout<<"You must re-root the phylogeny first"<<endl;
    throw 0;
  }
  BitSet gapSymbols=alignment.getGapSymbols();
  if(Q.getMatrixType()==MT_GAP) gapSymbols.purge();
  LCO_Likelihooder L(largestID+1,numAlpha,alphabet,alignment,column,
		     root,order,Q.getAlphabetMap(),gapSymbols);
  phylogeny.postorderTraversal(L);
  return L.getLikelihood();
}



// -----------------------------------------------------------------
// log-likelihood summed over all columns
// -----------------------------------------------------------------
double LCO_Felsenstein::logLikelihood()
{
  if(numThreads>1)
    return logLikelihood_threaded();
  
  double L=0;
  int numColumns=alignment.getLength();
  for(int i=0 ; i<numColumns ; ++i)
    L+=logLikelihood(i);
  return L;
}



// -----------------------------------------------------------------
// log-likelihood for a range of columns
// -----------------------------------------------------------------
double LCO_Felsenstein::logLikelihood(int begin,int end)
{
  double L=0;
  for(int i=begin ; i<end ; ++i)
    L+=logLikelihood(i);
  return L;
}



double LCO_Felsenstein::logLikelihood3(int begin,int end)
{
  return logLikelihoodInPhase(begin,end,2);
}



double LCO_Felsenstein::logLikelihoodInPhase(int begin,int end,int phase)
{
  double L=0;
  for(int i=begin ; i<end ; ++i)
    if(i%3==phase)
      L+=logLikelihood(i);
  return L;
}



double LCO_Felsenstein::logLikelihood_bestFrame(int begin,int end)
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
double LCO_Felsenstein::logLikelihood_threaded()
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
        LCO_FelThread *thread=new LCO_FelThread(*this,begin,end);
        threads[i]=thread;
        thread->start();
        nextCol=end;
    }

    // Wait for all the threads to finish
    double L=0;
    for(int i=0 ; i<numThreads ; ++i)
    {
        LCO_FelThread *thread=threads[i];
        thread->join();
        L+=thread->getLogLikelihood();
        delete thread;
    }

    // Return the likelihood
    return L;
}



/****************************************************************
                           LCO_FelThread methods
 ****************************************************************/
LCO_FelThread::LCO_FelThread(LCO_Felsenstein &alg,int begin,int end)
    : alg(alg),
      begin(begin),
      end(end)
{
    // ctor
}



void LCO_FelThread::f()
{
    L=0.0;
    for(int i=begin ; i<end ; ++i)
        L+=alg.logLikelihood(i);
}



double LCO_FelThread::getLogLikelihood()
{
    return L;
}



