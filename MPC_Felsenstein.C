/****************************************************************
 MPC.C
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
#include "MPC.H"
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
    Symbol N;
    bool skipN;
    RootNode *root; // ### for debugging only
    
    FitchLikelihooder(int numNodes,int numAlpha,Alphabet &alpha,
                      MultSeqAlignment &A,int column,bool skipN,
                      RootNode *root,int order)
        : A(A), numNodes(numNodes), numAlpha(numAlpha),
          alpha(alpha), L(numNodes,numAlpha), column(column),
          N(alpha.lookup('N')), skipN(skipN), order(order),
          root(root), likelihood(NEGATIVE_INFINITY)
        {
	  L.setAllTo(NEGATIVE_INFINITY);
        }
    
    virtual void processNode(LeafNode &u) // DO NOT CHANGE THIS FUNCTION
        {
            int id=u.getID();
            Symbol a=A.getIthTrack(id)[column];
            Array2D<double>::RowIn2DArray<double> row=L[id];
            if(a==N)
                for(Symbol i=0 ; i<numAlpha; ++i)
                    row[i]=log(1);// missing data -- same as Seipel & Haussler
            else
                for(Symbol i=0 ; i<numAlpha; ++i) 
                    row[i]=(i==a ? log(1) : log(0));
        }
    
    virtual void processNode(RootNode &u) {
        int id=u.getID();
        Symbol a=A.getIthTrack(id)[column];
        int childID=u.getChild()->getID();
        NthOrdSubstMatrix &Pt=*u.getSubstMatrix();
        if(a==N) {
	  likelihood=0; return;// ### ignore gaps in target sequence

	  const double prior=0.25; // ### need to replace this with true prior
	  Vector<double> V;
	  for(Symbol a=0 ; a<numAlpha ; ++a)
	    if(a!=N) V.push_back(prior*processInternalChild(a,id,childID,Pt));
	  likelihood=sumLogProbs(V);
        }
        else
            likelihood=processInternalChild(a,id,childID,Pt);
    }
    
    virtual void processNode(InternalNode &u) {
        int id=u.getID();
	NthOrdSubstMatrix &leftPt=*u.getLeftSubstMatrix();
	NthOrdSubstMatrix &rightPt=*u.getRightSubstMatrix();
	int left=u.getLeft()->getID(), right=u.getRight()->getID();
	Array2D<double>::RowIn2DArray<double> row=L[id];
        Symbol a=A.getIthTrack(id)[column];
	if(a==N)
	  for(a=0 ; a<numAlpha ; ++a) {
            if(skipN && a==N) continue;
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
        track.getSeq().getSubsequence(contextBegin,contextLen,context);
	int pos=MultSeqAlignment::rightmostGapPos(context,N);
	if(pos>=0) {
	  Sequence temp;
	  context.getSubsequence(pos+1,contextLen-pos-1,temp);
	  context=temp;
	}
        SubstitutionMatrix &Pt=*noPt.lookup(context,0,-1);
        Array2D<double>::RowIn2DArray<double> row=L[child];
        Symbol childSym=childTrack.getSeq()[column];
        double ll;
	if(childSym==N){
	  Vector<double> V;
	  for(Symbol a=0 ; a<numAlpha ; ++a) {
            if(a!=N) {
	      V.push_back(safeAdd(row[a],log(Pt(parentSymbol,a))));
	    }
	  }
	  ll=sumLogProbs(V);
	}
	else
	  ll=safeAdd(row[childSym],log(Pt(parentSymbol,childSym)));
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
MPC::MPC(int order,
				     Phylogeny &P,
				     NthOrdRateMatrix &Q,
				     MultSeqAlignment &A,
				     MolecularSequenceType T,
				     int numThreads)
    : phylogeny(P),
      Q(Q),
      alignment(A),
      seqType(T),
      alphabet(T==DNA ? (Alphabet&) DnaAlphabet::global : 
               (Alphabet&) AminoAlphabet::global),
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
void MPC::attachPhylogenyNodes()
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
double MPC::logLikelihood(int column)
{
  PhylogenyNode *root=phylogeny.getRoot();
  if(root->getNodeType()!=ROOT_NODE) {
    cout<<"You must re-root the phylogeny first"<<endl;
    throw 0;
  }
  //int rootID=root->getID();
  FitchLikelihooder L(largestID+1,numAlpha,alphabet,alignment,column,
		 seqType==DNA,root,order);
  phylogeny.postorderTraversal(L);
  return L.getLikelihood();
}



// -----------------------------------------------------------------
// log-likelihood summed over all columns
// -----------------------------------------------------------------
double MPC::logLikelihood()
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
double MPC::logLikelihood(int begin,int end)
{
  double L=0;
  for(int i=begin ; i<end ; ++i)
    L+=logLikelihood(i);
  return L;
}



// -----------------------------------------------------------------
// multi-threaded version
// -----------------------------------------------------------------
double MPC::logLikelihood_threaded()
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
FitchThread::FitchThread(MPC &alg,int begin,int end)
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



