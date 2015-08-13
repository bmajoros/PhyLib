/****************************************************************
 SingleGainParsimony.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_SingleGainParsimony_H
#define INCL_SingleGainParsimony_H
#include <iostream>
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Array1D.H"
#include "BOOM/BitSet.H"
#include "Phylogeny.H"
using namespace std;
using namespace BOOM;


struct DolloUp;
struct DolloDown;


/****************************************************************
 Performs gap-pattern reconstruction over a phylogeny using the 
 Single-Gain Parsimony algorithm, which finds the most parsimonious
 gap/residue assignment within an alignment column, subject to the
 constraint that the column can have at most one insertion event.
 This is in keeping with the definition of orthology --- multiple
 characters arising independently (even if appearing identical and/
 or appearing at seemingly orthologous sites) cannot be considered
 orthologous, since they do not share a common ancestor.

 This is known in the literature as DOLLO PARSIMONY.

 ****************************************************************/
class SingleGainParsimony
{
public:
  SingleGainParsimony(Phylogeny &,const Alphabet &,MultSeqAlignment &,
		       const BitSet &gapSymbols);
  void run(int column);
  void run(); // all columns
private:
  Phylogeny &phylogeny;
  int numTaxa;
  const Alphabet &alphabet;
  int alphabetSize;
  MultSeqAlignment &alignment;
  Array1D<BitSet> sets;
  Array1D<AlignmentSeq*> tracks;
  BitSet gapSymbols;

  void up(int column);
  void down(int column);

  friend class DolloUp;
  friend class DolloDown;
};



/****************************************************************
            Utility classes used by SingleGainParsimony
 ****************************************************************/

class DolloUp : public TreeVisitor {
public:
  DolloUp(SingleGainParsimony &,int column,const BitSet &gapSymbols);
  virtual void processNode(InternalNode &);
  virtual void processNode(LeafNode &);
  virtual void processNode(RootNode &);
private:
  SingleGainParsimony &master;
  int column;
  BitSet gapSymbols;
};


class DolloDown : public TreeVisitor {
public:
  DolloDown(SingleGainParsimony &,int column,const BitSet &gapSymbols);
  virtual void processNode(InternalNode &);
  virtual void processNode(LeafNode &);
  virtual void processNode(RootNode &);
private:
  SingleGainParsimony &master;
  int column;
  BitSet gapSymbols;
  Symbol gapSymbol;
};


#endif

