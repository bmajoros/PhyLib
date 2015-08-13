/****************************************************************
 ACO_Felsenstein.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "ACO_Felsenstein.H"
#include "BOOM/DnaAlphabet.H"
#include "NthOrdSubstMatrix.H"
#include "BOOM/Array2D.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/BitSet.H"
#include "SafeAdd.H"
using namespace std;
using namespace BOOM;


/*
     NOTE:

     AS CURRENTLY IMPLEMENTED, AN Nth ORDER MODEL CAN'T PROCESS A
     SEQUENCE OF LENGTH SHORTER THAN N-1.  TO FIX THIS, MODIFY THE
     PROCESSING OF THE ROOT TAXON TO DETECT WHEN THE WINDOW EXTENDS
     PAST THE END OF THE SEQUENCE AND THEN ITERATE OVER ALL POSSIBLE
     EXTENSIONS OF THE SEQUENCE TO FILL THE WINDOW.
 */



// -----------------------------------------------------------------
// class ACO_Likelihooder
// -----------------------------------------------------------------

struct ACO_Likelihooder : public HOG_Likelihooder
{
  int order;
  NmerRateMatrix *Q;
  BitSet gapSymbols;
  Array1D<double> V;
  bool truncate;

public:
  ACO_Likelihooder(int numNodes,int numAlpha,Alphabet &alpha,
		    MultSeqAlignment &A,int firstCol,int lastCol,
		    RootNode *root,AlphabetMap &alphabetMap,
		    MolecularSequenceType seqType,int order,
		    NmerRateMatrix *Q,AlignmentNmerTable *nmerTable,
		    bool truncate)
    : HOG_Likelihooder(numNodes,numAlpha,alpha,A,firstCol,lastCol,
		       root,alphabetMap,seqType,nmerTable),
      order(order),
      Q(Q),
      gapSymbols(A.getGapSymbols()),
      truncate(truncate)
  {
    // ctor
    numNmers=pow((float)alphabetMap.getRangeSize(),(float)order+1);
    L.resize(numNodes,numNmers);
    V.resize(numNmers);
  }

  virtual void processNode(LeafNode &u) {
    int id=u.getID();
    Array2D<double>::RowIn2DArray<double> row=L[id];
    /*
    Sequence &track=A.getIthTrack(id).getSeq();
    Sequence leafNmer, nmer;
    track.getSubsequence(firstCol,numCols,leafNmer);
    for(int i=0 ; i<numNmers ; ++i) {
      nmer.fromInt(i,order+1,alphabetMap);
      row[i]=degenerateDnaMatch(nmer,leafNmer,gapSymbols,numCols) ? 
        0 : NEGATIVE_INFINITY;
    }
    */
    row.setAllTo(NEGATIVE_INFINITY);
    NmerProfile *profile=nmerTable->getNmer(id,firstCol,order);
    Vector<int>::iterator cur, end;
    if(truncate) {
      cur=profile->degenerateTruncated.begin();
      end=profile->degenerateTruncated.end();
    }
    else {
      cur=profile->degenerateMatches.begin();
      end=profile->degenerateMatches.end();
    }
    for(; cur!=end ; ++cur) row[*cur]=0;
  }
    
  virtual void processNode(RootNode &u) {
    // PRECONDITION: the root track contains no gap symbols!
    int id=u.getID();
    int childID=u.getChild()->getID();
    NthOrdSubstMatrix &Pt=*u.getSubstMatrix();
    /*
    Sequence &track=A.getIthTrack(id).getSeq();
    Sequence rootNmer;
    track.getSubsequence(firstCol,order+1,rootNmer);
    int rootNmerCode=rootNmer.asInt(alphabetMap);
    */
    NmerProfile *profile=nmerTable->getNmer(id,firstCol,order);
    int rootNmerCode=profile->ungappedNmerCode;
    likelihood=processInternalChild(rootNmerCode,childID,Pt);
  }
    
  virtual void processNode(InternalNode &u) {
    int id=u.getID();
    Array2D<double>::RowIn2DArray<double> row=L[id];
    int left=u.getLeft()->getID(), right=u.getRight()->getID();
    NthOrdSubstMatrix &leftPt=*u.getLeftSubstMatrix();
    NthOrdSubstMatrix &rightPt=*u.getRightSubstMatrix();
    for(int i=0 ; i<numNmers ; ++i) {
      row[i]=
	processInternalChild(i,left,leftPt)+
	processInternalChild(i,right,rightPt);
    }
  }

  virtual double processInternalChild(int parentNmerCode,
				      int childID,
				      NthOrdSubstMatrix &dualPt) {
    NmerSubstMatrix &Pt=static_cast<NmerSubstMatrix&>(dualPt);
    Array2D<double>::RowIn2DArray<double> row=L[childID];
    for(int i=0 ; i<numNmers ; ++i) {
      double joint=Pt.lookupNmers(parentNmerCode,i);
      V[i]=safeAdd(joint,row[i]);
    }
    double ll=sumLogProbs(V);
    return ll;
  }
};



// -----------------------------------------------------------------
// class ACO_Felsenstein
// -----------------------------------------------------------------

ACO_Felsenstein::ACO_Felsenstein(Phylogeny &phylogeny,
				 NmerRateMatrix *Q,
				 MultSeqAlignment &A,
				 AlignmentNmerTable *nmerTable,
				 int numThreads) 
  : HOG_Felsenstein(Q->getOrder(),phylogeny,*Q,A,DNA,nmerTable,numThreads)
{
  // ctor
}



HOG_Likelihooder *ACO_Felsenstein::getLikelihooder(int numNodes,int numAlpha,
		    Alphabet &alpha,MultSeqAlignment &A,int firstCol,
		    int lastCol,RootNode *root,
		    AlphabetMap &alphabetMap,MolecularSequenceType seqType)
{
  NmerRateMatrix *nmerQ=static_cast<NmerRateMatrix*>(&Q);
  bool truncate=lastCol-firstCol<order;
  return new ACO_Likelihooder(numNodes,numAlpha,alpha,A,firstCol,lastCol,
			       root,alphabetMap,seqType,order,nmerQ,nmerTable,
			       truncate);
}



