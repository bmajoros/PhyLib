/****************************************************************
 MUD.C
 Bill Majoros - bmajoros@duke.edu
 ****************************************************************/
#include "MUD.H"
#include <iostream>
#include "tigr++/TigrStack.H"
using namespace std;


MUD::MUD(Phylogeny &phylogeny)
  : phylogeny(phylogeny),
    leafNodes(phylogeny.getRoot()->getNumLeaves()),
    internalNodes(phylogeny.getRoot()->getNumLeaves()-1),
    numLeafNodes(0),
    numInternalNodes(0)
{
  // ctor

  root=(RootNode*) phylogeny.getRoot();
  preprocess(root);
}



void MUD::preprocess(PhylogenyNode *node)
{
  switch(node->getNodeType())
    {
    case ROOT_NODE:
      {
	RootNode *root=(RootNode*) node;
	preprocess(root->getChild());
      }
      break;
    case LEAF_NODE:
      leafNodes[numLeafNodes++]=(LeafNode*) node;
      break;
    case INTERNAL_NODE:
      {
	InternalNode *internalNode=(InternalNode*) node;
	preprocess(internalNode->getLeft());
	preprocess(internalNode->getRight());
	internalNodes[numInternalNodes++]=internalNode;
      }
      break;
    }
}



void MUD::go(MultiAlignment &alignment)
{
  char c;

  // Process the entire length of the alignment
  int alignmentLength=alignment.getLength();
  for(int pos=0 ; pos<alignmentLength ; ++pos)
    {
      // Process all leaf nodes for this position
      for(int i=0 ; i<numLeafNodes ; ++i)
	{
	  LeafNode *node=leafNodes[i];
	  ResidueCounter &counter=node->getCounter();
	  counter.clear();
	  c=node->getResidue(pos);
	  //if(c==' ' || c=='-') continue;
	  if(c==' ') continue;
	  counter.increment(c);
	}

      // Process all internal nodes for this position
      for(int i=0 ; i<numInternalNodes ; ++i) 
	internalNodes[i]->clearResolvedBit();
      for(int i=0 ; i<numInternalNodes ; ++i)
	{
	  InternalNode *node=internalNodes[i];
	  ResidueCounter &thisCounter=node->getCounter();
	  thisCounter.clear();
	  ResidueCounter &leftCounter=node->getLeft()->getCounter();
	  ResidueCounter &rightCounter=node->getRight()->getCounter();
	  leftCounter.add(rightCounter,thisCounter);
	  if(node->getParent()==root)
	    {
	      c=(*root->getAlignmentTrack())[pos];
	      thisCounter.increment(c);
	    }
	  if(thisCounter.getMostPopular(c) ||
	     node->getParent()==root)
	    {
	      (*node->getAlignmentTrack())[pos]=c;
	      node->setResolvedBit();
	      thisCounter.clear();
	      thisCounter.increment(c);
	      resolve(node->getLeft(),pos,c);
	      resolve(node->getRight(),pos,c);
	    }
	}
    }
}



void MUD::resolve(PhylogenyNode *node,int pos,char c)
{
  if(node->getNodeType()!=INTERNAL_NODE) return;
  InternalNode *iNode=(InternalNode*) node;
  if(iNode->isResolved()) return;
  (*iNode->getAlignmentTrack())[pos]=c;
  resolve(iNode->getLeft(),pos,c);
  resolve(iNode->getRight(),pos,c);
}
