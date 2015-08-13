/****************************************************************
 RCO_Felsenstein.C
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
#include "RCO_Felsenstein.H"
using namespace std;
using namespace BOOM;


// -----------------------------------------------------------------
// This is where the actual likelihood algorithm is implemented...
// -----------------------------------------------------------------
struct RCO_Likelihooder : public TreeVisitor
{
public:
  Array2D<double> L; // the dynamic programming matrix
  MultSeqAlignment &A;
  int numNodes, numAlpha, column;
  Alphabet &alpha;
  AlphabetMap &alphabetMap;
  double likelihood;
  Sequence rootContext, noContext;
  int rootContextCode, rootContextLength;
  BitSet gapSymbols;

  RCO_Likelihooder(int numNodes,int numAlpha,Alphabet &alpha,
		   MultSeqAlignment &A,int column,
		   RootNode *root,int order,AlphabetMap &alphabetMap,
		   const BitSet &gapSymbols)
    : A(A), numNodes(numNodes), numAlpha(numAlpha),
      alpha(alpha), L(numNodes,numAlpha), column(column),
      alphabetMap(alphabetMap), gapSymbols(gapSymbols)
  {
    // ctor
    
    // Get the root's left-context for this column (for RCO)
    int rootID=root->getID();
    AlignmentSeq &rootTrack=A.getIthTrack(rootID);
    int contextBegin=column-order;
    if(contextBegin<0) contextBegin=0;
    int contextLen=column-contextBegin;
    rootTrack.getSeq().getSubsequence(contextBegin,contextLen,
				      rootContext);
    int pos=MultSeqAlignment::rightmostGapPos(rootContext,gapSymbols);
    if(pos>=0) {
      Sequence temp;
      rootContext.getSubsequence(pos+1,contextLen-pos-1,temp);
      rootContext=temp;
    }
    rootContextCode=rootContext.asInt(alphabetMap);
    rootContextLength=contextLen;
  }
  
  virtual void processNode(LeafNode &u) // DO NOT CHANGE THIS FUNCTION
  {
    int id=u.getID();
    Symbol a=A.getIthTrack(id)[column];
    Array2D<double>::RowIn2DArray<double> row=L[id];
    if(gapSymbols.isMember(a))
      for(Symbol i=0 ; i<numAlpha; ++i)
	row[i]=0; //=log(1); missing data -- same as Seipel & Haussler
    else
      for(Symbol i=0 ; i<numAlpha; ++i) 
	row[i]=(i==a ? 0 : NEGATIVE_INFINITY);
  }
  
  virtual void processNode(RootNode &u) {
    int id=u.getID();
    Symbol a=A.getIthTrack(id)[column];
    int childID=u.getChild()->getID();
    NthOrdSubstMatrix &Pt=*u.getSubstMatrix();
    if(gapSymbols.isMember(a))
      likelihood=0; // shouldn't happen...
    else
      likelihood=processInternalChild(a,childID,Pt,rootContextCode,
				      rootContextLength);
  }
  
  virtual void processNode(InternalNode &u) {
    int id=u.getID();
    int left=u.getLeft()->getID(), right=u.getRight()->getID();
    NthOrdSubstMatrix &leftPt=*u.getLeftSubstMatrix();
    NthOrdSubstMatrix &rightPt=*u.getRightSubstMatrix();
    Array2D<double>::RowIn2DArray<double> row=L[id];
    for(Symbol a=0 ; a<numAlpha ; ++a) {
      if(gapSymbols.isMember(a)) continue;
      row[a]=
	processInternalChild(a,left,leftPt,0,0)+
	processInternalChild(a,right,rightPt,0,0);
    }
  }
  
  double processInternalChild(Symbol parentSymbol,int child,
			      NthOrdSubstMatrix &noPt,
			      int contextCode,int contextLength) {
    SubstitutionMatrix &Pt=*noPt.lookup(contextCode,contextLength);
    Array2D<double>::RowIn2DArray<double> row=L[child];
    Vector<double> V;
    for(Symbol b=0 ; b<numAlpha ; ++b)
      if(!gapSymbols.isMember(b))
	V.push_back(row[b]+Pt(parentSymbol,b)); // Pt already in log space
    double ll=sumLogProbs(V);
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
RCO_Felsenstein::RCO_Felsenstein(int order,
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
void RCO_Felsenstein::attachPhylogenyNodes()
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
double RCO_Felsenstein::logLikelihood(int column)
{
  PhylogenyNode *root=phylogeny.getRoot();
  if(root->getNodeType()!=ROOT_NODE) {
    cout<<"You must re-root the phylogeny first"<<endl;
    throw 0;
  }
  RCO_Likelihooder L(largestID+1,numAlpha,alphabet,alignment,column,
		     root,order,Q.getAlphabetMap(),alignment.getGapSymbols());
  phylogeny.postorderTraversal(L);
  return L.getLikelihood();
}



// -----------------------------------------------------------------
// log-likelihood summed over all columns
// -----------------------------------------------------------------
double RCO_Felsenstein::logLikelihood()
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
double RCO_Felsenstein::logLikelihood(int begin,int end)
{
  double L=0;
  for(int i=begin ; i<end ; ++i)
    L+=logLikelihood(i);
  return L;
}


double RCO_Felsenstein::logLikelihood3(int begin,int end)
{
  return logLikelihoodInPhase(begin,end,2);
}



double RCO_Felsenstein::logLikelihoodInPhase(int begin,int end,int phase)
{
  double L=0;
  for(int i=begin ; i<end ; ++i)
    if(i%3==phase)
      L+=logLikelihood(i);
  return L;
}



double RCO_Felsenstein::logLikelihood_bestFrame(int begin,int end)
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
double RCO_Felsenstein::logLikelihood_threaded()
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
        RCO_FelThread *thread=new RCO_FelThread(*this,begin,end);
        threads[i]=thread;
        thread->start();
        nextCol=end;
    }

    // Wait for all the threads to finish
    double L=0;
    for(int i=0 ; i<numThreads ; ++i)
    {
        RCO_FelThread *thread=threads[i];
        thread->join();
        L+=thread->getLogLikelihood();
        delete thread;
    }

    // Return the likelihood
    return L;
}



/****************************************************************
                           RCO_FelThread methods
 ****************************************************************/
RCO_FelThread::RCO_FelThread(RCO_Felsenstein &alg,int begin,int end)
    : alg(alg),
      begin(begin),
      end(end)
{
    // ctor
}



void RCO_FelThread::f()
{
    L=0.0;
    for(int i=begin ; i<end ; ++i)
        L+=alg.logLikelihood(i);
}



double RCO_FelThread::getLogLikelihood()
{
    return L;
}



