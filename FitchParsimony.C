/****************************************************************
 FitchParsimony.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "FitchParsimony.H"
#include "BOOM/Random.H"
using namespace std;
using namespace BOOM;

/*
void printSet(BitSet &bs,const Alphabet &alpha) {
    int n=bs.cardinality();
    cout<<"{";
    for(int i=0 ; i<n ; ++i) {
        Symbol s=bs.getIthMember(i);
        cout<<alpha.lookup(s);
        if(i<n-1) cout<<",";
    }
    cout<<"}";
}
*/


FitchParsimony::FitchParsimony(Phylogeny &phylogeny,const Alphabet &alphabet,
                               MultSeqAlignment &alignment,
			       const BitSet &gapSymbols)
    : phylogeny(phylogeny),
      alphabet(alphabet),
      alignment(alignment),
      numTaxa(phylogeny.getNumNodes()),
      alphabetSize(alphabet.size()),
      gapSymbols(gapSymbols)
{
    // ctor

    initialSets.resize(numTaxa);
    finalSets.resize(numTaxa);
    tracks.resize(numTaxa);
    for(int i=0 ; i<numTaxa ; ++i) {
        initialSets[i].setSize(alphabetSize);
        finalSets[i].setSize(alphabetSize);
    }
    TrackCollector trackCollector(alignment,tracks);
    phylogeny.postorderTraversal(trackCollector);
}



void FitchParsimony::run() {
  int numColumns=alignment.getLength();
  for(int i=0 ; i<numColumns ; ++i)
    run(i);
}



void FitchParsimony::run(int column) {
  for(int i=0 ; i<numTaxa ; ++i) {
    initialSets[i].purge();
    finalSets[i].purge();
  }
  up(column);
  down(column);
}



void FitchParsimony::up(int column) {
  FitchUp fitchUp(*this,column,gapSymbols);
  phylogeny.postorderTraversal(fitchUp);
}



void FitchParsimony::down(int column) {
  FitchDown fitchDown(*this,column,gapSymbols);
  phylogeny.preorderTraversal(fitchDown);
}


/****************************************************************
                          FitchUp methods
 ****************************************************************/
FitchUp::FitchUp(FitchParsimony &master,int column,const BitSet &gapSymbols)
    : master(master),
      column(column),
      gapSymbols(gapSymbols)
{
    // ctor
}


void FitchUp::processNode(InternalNode &V) {
  PhylogenyNode *leftNode=V.getLeft(), *rightNode=V.getRight();
  int ID=V.getID(), leftID=leftNode->getID(), rightID=rightNode->getID();
  BitSet &leftSet=master.initialSets[leftID];
  BitSet &rightSet=master.initialSets[rightID];
  BitSet &thisSet=master.initialSets[ID];
  leftSet.intersect(rightSet,thisSet);
  if(thisSet.cardinality()==0) leftSet.unionWith(rightSet,thisSet);
  //if(V.getParent()==NULL) printSet(thisSet,master.alphabet);
}



void FitchUp::processNode(LeafNode &V) {
    int id=V.getID();
    AlignmentSeq &track=*master.tracks[id];
    Symbol s=track.getSeq()[column];
    BitSet &bitSet=master.initialSets[id];
    if(!gapSymbols.isMember(s)) bitSet.addMember(s);
}



void FitchUp::processNode(RootNode &V) {
    int id=V.getID();
    AlignmentSeq &track=*master.tracks[id];
    Symbol s=track.getSeq()[column];
    master.initialSets[id].addMember(s);
    master.finalSets[id].addMember(s);
}



/****************************************************************
                         FitchDown methods
 ****************************************************************/
FitchDown::FitchDown(FitchParsimony &master,int column,
		     const BitSet &gapSymbols)
    : master(master),
      column(column),
      gapSymbols(gapSymbols),
      gapSymbol(INVALID_SYMBOL)
{
  int n=gapSymbols.getMaxSize();
  for(int i=0 ; i<n ; ++i)
    if(gapSymbols.isMember(i))
      gapSymbol=i;
}


void FitchDown::processNode(InternalNode &V) {
  int ID=V.getID();
  BitSet &Fp=master.finalSets[ID];
  BitSet &Sp=master.initialSets[ID];
  PhylogenyNode *parent=V.getParent();
  if(parent) {
    PhylogenyNode *left=V.getLeft(), *right=V.getRight();
    int parentID=parent->getID(), leftID=left->getID(), rightID=right->getID();
    BitSet &Fa=master.finalSets[parentID];
    Sp.intersect(Fa,Fp);
    if(Fp!=Fa) {
      BitSet &Sq=master.initialSets[leftID];
      BitSet &Sr=master.initialSets[rightID];
      BitSet Sq_I_Sr(master.alphabetSize);
      Sq.intersect(Sr,Sq_I_Sr);
      if(Sq_I_Sr.cardinality()>0) {
	BitSet Sq_U_Sr(master.alphabetSize);
	Sq.unionWith(Sr,Sq_U_Sr);
	BitSet Sq_U_Sr_I_Fa(master.alphabetSize);
	Sq_U_Sr.intersect(Fa,Sq_U_Sr_I_Fa);
	Sq_U_Sr_I_Fa.unionWith(Sp,Fp);
      }
      else {
	Sp.unionWith(Fa,Fp);
      }
    }
  }
  else Fp=Sp;
  
  // Resolve ties and store the inferred symbols into the alignment
  int n=Fp.cardinality();
  if(n>1) {
    int element=0;//###RandomNumber(n); -> randomness is bad for training!!!
    Symbol s=Fp.getIthMember(element);
    Fp.purge();
    Fp.addMember(s);
  }
  Symbol s=n>0 ? Symbol(Fp.getIthMember(0)) : gapSymbol;

  (*master.tracks[ID])[column]=s;
}



void FitchDown::processNode(LeafNode &V) {
    // do nothing
}



void FitchDown::processNode(RootNode &V) {
    // do nothing
}



/****************************************************************
                      TrackCollector methods
 ****************************************************************/
TrackCollector::TrackCollector(MultSeqAlignment &alignment,
                               Array1D<AlignmentSeq*> &tracks)
    : alignment(alignment),
      tracks(tracks)
{
    // ctor
}



void TrackCollector::processNode(InternalNode &V) {process(V);}



void TrackCollector::processNode(LeafNode &V) {process(V);}



void TrackCollector::processNode(RootNode &V) {process(V);}



void TrackCollector::process(PhylogenyNode &V) {
    int id=V.getID();
    const String &name=V.getName();
    AlignmentSeq *track=&alignment.findOrCreateTrack(name);
    int L=alignment.getLength();
    track->extendToLength(L,alignment.getGapSymbol());
    tracks[id]=track;
}





