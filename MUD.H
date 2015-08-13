/****************************************************************
 MUD.H : the Majority Up-Down algorithm for ancestor inference
 Bill Majoros - bmajoros@duke.edu
 ****************************************************************/
#ifndef INCL_MUD_H
#define INCL_MUD_H
#include "tigr++/TigrArray1D.H"
#include "tigr++/MultiAlignment.H"
#include "Phylogeny.H"

class MUD
{
  TigrArray1D<LeafNode*> leafNodes;
  TigrArray1D<InternalNode*> internalNodes;
  int numLeafNodes, numInternalNodes;
  Phylogeny &phylogeny;
  RootNode *root;

  void preprocess(PhylogenyNode*);
  inline void resolve(PhylogenyNode *,int pos,char c);
public:
  MUD(Phylogeny &);
  void go(MultiAlignment &);
};


#endif
