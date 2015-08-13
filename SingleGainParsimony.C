/****************************************************************
 SingleGainParsimony.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "SingleGainParsimony.H"
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


SingleGainParsimony::SingleGainParsimony(Phylogeny &phylogeny,
					   const Alphabet &alphabet,
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
  
  sets.resize(numTaxa);
  tracks.resize(numTaxa);
  for(int i=0 ; i<numTaxa ; ++i)
    sets[i].setSize(alphabetSize);
  TrackCollector trackCollector(alignment,tracks);
  phylogeny.postorderTraversal(trackCollector);
}



void SingleGainParsimony::run() {
  int numColumns=alignment.getLength();
  for(int i=0 ; i<numColumns ; ++i)
    run(i);
}



void SingleGainParsimony::run(int column) {
  for(int i=0 ; i<numTaxa ; ++i)
    sets[i].purge();
  up(column);
  down(column);
}



void SingleGainParsimony::up(int column) {
  DolloUp dolloUp(*this,column,gapSymbols);
  phylogeny.postorderTraversal(dolloUp);
}



void SingleGainParsimony::down(int column) {
  DolloDown dolloDown(*this,column,gapSymbols);
  phylogeny.preorderTraversal(dolloDown);
}


/****************************************************************
                          DolloUp methods
 ****************************************************************/
DolloUp::DolloUp(SingleGainParsimony &master,int column,const BitSet 
	       &gapSymbols)
  : master(master), column(column), gapSymbols(gapSymbols)
{
  // ctor
}



void DolloUp::processNode(InternalNode &V) {
  PhylogenyNode *leftNode=V.getLeft(), *rightNode=V.getRight();
  int ID=V.getID(), leftID=leftNode->getID(), rightID=rightNode->getID();
  BitSet &leftSet=master.sets[leftID], &rightSet=master.sets[rightID];
  BitSet &thisSet=master.sets[ID];
  leftSet.unionWith(rightSet,thisSet);
}



void DolloUp::processNode(LeafNode &V) {
  int id=V.getID();
  AlignmentSeq &track=*master.tracks[id];
  Symbol s=track.getSeq()[column];
  if(!gapSymbols.isMember(s)) master.sets[id].addMember(s);
}



void DolloUp::processNode(RootNode &V) {
  // do nothing
}



/****************************************************************
                         DolloDown methods
 ****************************************************************/
DolloDown::DolloDown(SingleGainParsimony &master,int column,
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


void DolloDown::processNode(InternalNode &V) {
  int nAlpha=master.alphabetSize;
  PhylogenyNode *left=V.getLeft(), *right=V.getRight();
  int ID=V.getID(), leftID=left->getID(), rightID=right->getID();
  BitSet &L=master.sets[leftID], &R=master.sets[rightID];
  BitSet &X=master.sets[ID];
  BitSet LiR(nAlpha), PiX(nAlpha);
  L.intersect(R,LiR);
  PhylogenyNode *parent=V.getParent();
  if(parent) {
    master.sets[parent->getID()].intersect(X,PiX);
    LiR.unionWith(PiX,X);
  }
  else X=LiR;
  
  // Resolve ties and store the inferred symbols into the alignment
  int n=X.cardinality();
  if(n>1) {
    Symbol s=X.getIthMember(0); // don't pick randomly -- bad for derivatives
    X.purge();
    X.addMember(s);
  }
  Symbol s=n>0 ? Symbol(X.getIthMember(0)) : gapSymbol;
  //if(s==Symbol(0)) throw "GPE_UNKNOWN IN DOLLO";

  (*master.tracks[ID])[column]=s;
}



void DolloDown::processNode(LeafNode &V) {
    // do nothing
}



void DolloDown::processNode(RootNode &V) {
    // do nothing
}

