/****************************************************************
 FelsensteinsAlgorithm.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/DnaAlphabet.H"
#include "BOOM/AminoAlphabet.H"
#include "BOOM/Array2D.H"
#include "FelsensteinsAlgorithm.H"
using namespace std;
using namespace BOOM;


inline double safeAdd(double a,double b)
{
    if(isinf(a) || isinf(b)) return a;
    else return a+b;
}



struct Likelihooder : public TreeVisitor
{
public:
    Array2D<double> L;
    MultiAlignment &A;
    int numNodes, numAlpha, column;
    Alphabet &alpha;
    double likelihood;
    Symbol N;
    bool skipN;
    
    Likelihooder(int numNodes,int numAlpha,Alphabet &alpha,
		 MultiAlignment &A,int column,bool skipN)
        : A(A), numNodes(numNodes), numAlpha(numAlpha),
          alpha(alpha), L(numNodes,numAlpha), column(column),
          N(alpha.lookup('N')), skipN(skipN) {}
    virtual void processNode(LeafNode &u)
        {
            int id=u.getID();
            char c=A.getIthTrack(id)[column];
            if(!alpha.isDefined(c)) c='N';// dashes & other non-ACGT chars
            Symbol a=alpha.lookup(c);
            Array2D<double>::RowIn2DArray<double> row=L[id];
            if(a==N)
                for(Symbol i=0 ; i<numAlpha; ++i)
                    row[i]=log(1);// missing data -- same as Seipel & Haussler
            else
                for(Symbol i=0 ; i<numAlpha; ++i) 
                    row[i]=(i==a ? log(1) : log(0));
        }
    double processInternalChild(Symbol parentSymbol,int child,
				SubstitutionMatrix &Pt)
        {
            Array2D<double>::RowIn2DArray<double> row=L[child];
            Array1D<double> V(numAlpha);
            for(Symbol b=0 ; b<numAlpha ; ++b)
                if(b!=N || !skipN)
                    V[b]=row[b]+log(Pt(parentSymbol,b));
            Symbol gamma=0;
            for(Symbol b=1 ; b<numAlpha ; ++b)
                if(b!=N || !skipN)
                    if(V[b]>V[gamma]) gamma=b;
            double Vgamma=V[gamma];
            double sum=0;
            for(Symbol b=0 ; b<numAlpha ; ++b)
                if(b!=N || !skipN)
                    if(b!=gamma)
                        sum+=exp(V[b]-Vgamma);
            double ll=Vgamma+log(1+sum);
            return ll;
        }
    virtual void processNode(InternalNode &u) 
        {
            int id=u.getID();
            Array2D<double>::RowIn2DArray<double> row=L[id];
            for(Symbol a=0 ; a<numAlpha ; ++a)
            {
                if(skipN && a==N) continue;
                int left=u.getLeft()->getID(), right=u.getRight()->getID();
                SubstitutionMatrix &leftPt=*u.getLeftSubstMatrix();
                SubstitutionMatrix &rightPt=*u.getRightSubstMatrix();
                row[a]=
                    processInternalChild(a,left,leftPt)+
                    processInternalChild(a,right,rightPt);
            }
        }
    virtual void processNode(RootNode &u)
        {
            int id=u.getID();
            char c=A.getIthTrack(id)[column];
            if(!alpha.isDefined(c)) c='N';
            Symbol a=alpha.lookup(c);
            int childID=u.getChild()->getID();
            SubstitutionMatrix &Pt=*u.getSubstMatrix();
            if(c=='N')//### ?
            {
                likelihood=0;
                for(Symbol a=0 ; a<numAlpha ; ++a)
                    if(skipN && a==N) continue;
                    else likelihood+=processInternalChild(a,childID,Pt);
            }
            else
                likelihood=processInternalChild(a,childID,Pt);
        }
    double getLikelihood()
        {
            return likelihood;
        }
};



FelsensteinsAlgorithm::FelsensteinsAlgorithm(Phylogeny &P,
					     RateMatrix &Q,
					     MultiAlignment &A,
					     MolecularSequenceType T,
                                             int numThreads)
    : phylogeny(P),
      Q(Q),
      alignment(A),
      seqType(T),
      alphabet(T==DNA ? (Alphabet&) DnaAlphabet::global : 
               (Alphabet&) AminoAlphabet::global),
      numThreads(numThreads),
      threads(numThreads)
{
    numAlpha=alphabet.size();
    attachPhylogenyNodes();
}



void FelsensteinsAlgorithm::attachPhylogenyNodes()
{
    struct Attacher : public TreeVisitor
    {
        MultiAlignment &A;
        RateMatrix &Q;
        int numNodes, largestID;
        
        Attacher(MultiAlignment &A,RateMatrix &Q) 
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



double FelsensteinsAlgorithm::logLikelihood(int column)
{
    if(phylogeny.getRoot()->getNodeType()!=ROOT_NODE)
    {
        cout<<"You must re-root the phylogeny first"<<endl;
        throw 0;
    }
    Likelihooder L(largestID+1,numAlpha,alphabet,alignment,column,
                   seqType==DNA);
    phylogeny.postorderTraversal(L);
    return L.getLikelihood();
}



double FelsensteinsAlgorithm::logLikelihood()
{
  if(numThreads>1)
    return logLikelihood_threaded();
  
  double L=0;
  int numColumns=alignment.getLength();
  for(int i=0 ; i<numColumns ; ++i)
    L+=logLikelihood(i);
  return L;
}



double FelsensteinsAlgorithm::logLikelihood(int begin,int end)
{
  double L=0;
  for(int i=begin ; i<end ; ++i)
    L+=logLikelihood(i);
  return L;
}



double FelsensteinsAlgorithm::logLikelihood_threaded()
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
FelThread::FelThread(FelsensteinsAlgorithm &alg,int begin,int end)
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




/***********************************************************************




                              OBSOLETE
 



double FelsensteinsAlgorithm::likelihood(int column)
{
  struct Likelihooder : public TreeVisitor
  {
  public:
    Array2D<double> L;
    MultiAlignment &A;
    int numNodes, numAlpha, column;
    Alphabet &alpha;
    double likelihood;
    Symbol N;
    bool skipN;

    Likelihooder(int numNodes,int numAlpha,Alphabet &alpha,
		 MultiAlignment &A,int column,bool skipN)
      : A(A), numNodes(numNodes), numAlpha(numAlpha),
	alpha(alpha), L(numNodes,numAlpha), column(column),
	N(alpha.lookup('N')), skipN(skipN) {}
    virtual void processNode(LeafNode &u)
    {
      int id=u.getID();
      char c=A.getIthTrack(id)[column];
//      if(!alpha.isDefined(c)) c='N';//### including dashes...
      Symbol a=alpha.lookup(c);
      Array2D<double>::RowIn2DArray<double> row=L[id];
      for(Symbol i=0 ; i<numAlpha; ++i) 
          row[i]=(i==a ? 1 : 0);
    }
    double processInternalChild(Symbol a,int child,SubstitutionMatrix &Pt)
    {
      double sum=0;
      Array2D<double>::RowIn2DArray<double> row=L[child];
      for(Symbol b=0 ; b<numAlpha ; ++b)
	if(b!=N || !skipN)
	  sum+=row[b]*Pt(a,b);
      return sum;
    }
    virtual void processNode(InternalNode &u) 
    {
      int id=u.getID();
      Array2D<double>::RowIn2DArray<double> row=L[id];
      for(Symbol a=0 ; a<numAlpha ; ++a)
	{
	  if(skipN && a==N) continue;
	  int left=u.getLeft()->getID(), right=u.getRight()->getID();
	  SubstitutionMatrix &leftPt=*u.getLeftSubstMatrix();
	  SubstitutionMatrix &rightPt=*u.getRightSubstMatrix();
	  row[a]=
	    processInternalChild(a,left,leftPt)*
	    processInternalChild(a,right,rightPt);
	}
    }
    virtual void processNode(RootNode &u)
    {
      int id=u.getID();
      char c=A.getIthTrack(id)[column];
      Symbol a=alpha.lookup(c);
      int childID=u.getChild()->getID();
      likelihood=0;
      SubstitutionMatrix &Pt=*u.getSubstMatrix();
      for(Symbol b=0 ; b<numAlpha ; ++b)
	if(b!=N || !skipN)
	  likelihood+=L[childID][b]*Pt(a,b);
    }
    double getLikelihood()
    {
      return likelihood;
    }
  };

  Likelihooder L(largestID+1,numAlpha,alphabet,alignment,column,
		 seqType==DNA);
  phylogeny.postorderTraversal(L);
  return L.getLikelihood();
}
*/



/*
double FelsensteinsAlgorithm::likelihood()
{
  double L=0;
  int numColumns=alignment.getLength();
  for(int i=0 ; i<numColumns ; ++i)
    L+=log(likelihood(i));
  return exp(L);
}
*/


