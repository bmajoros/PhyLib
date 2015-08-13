/****************************************************************
 FitchParsimony.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_FitchParsimony_H
#define INCL_FitchParsimony_H
#include <iostream>
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Array1D.H"
#include "BOOM/BitSet.H"
#include "Phylogeny.H"
using namespace std;
using namespace BOOM;


struct FitchUp;
struct FitchDown;


/****************************************************************
 Implements Fitch's maximum parsimony algorithm for ancestral
 state reconstruction.
 ****************************************************************/
class FitchParsimony
{
public:
  FitchParsimony(Phylogeny &,const Alphabet &,MultSeqAlignment &,
		 const BitSet &gapSymbols);
    void run(int column);
    void run(); // all columns
private:
    Phylogeny &phylogeny;
    int numTaxa;
    const Alphabet &alphabet;
    int alphabetSize;
    MultSeqAlignment &alignment;
    Array1D<BitSet> initialSets, finalSets;
    Array1D<AlignmentSeq*> tracks;
    BitSet gapSymbols;

    void up(int column);
    void down(int column);

    friend class FitchUp;
    friend class FitchDown;
};



/****************************************************************
            Utility classes used by FitchParsimony
 ****************************************************************/

class TrackCollector : public TreeVisitor {
public:
    TrackCollector(MultSeqAlignment &,Array1D<AlignmentSeq*> &);
    virtual void processNode(InternalNode &);
    virtual void processNode(LeafNode &);
    virtual void processNode(RootNode &);
private:
    MultSeqAlignment &alignment;
    Array1D<AlignmentSeq*> &tracks;
    void process(PhylogenyNode &);
};


class FitchUp : public TreeVisitor {
public:
  FitchUp(FitchParsimony &,int column,const BitSet &gapSymbols);
  virtual void processNode(InternalNode &);
  virtual void processNode(LeafNode &);
  virtual void processNode(RootNode &);
private:
  FitchParsimony &master;
  int column;
  BitSet gapSymbols;
};


class FitchDown : public TreeVisitor {
public:
  FitchDown(FitchParsimony &,int column,const BitSet &gapSymbols);
  virtual void processNode(InternalNode &);
  virtual void processNode(LeafNode &);
  virtual void processNode(RootNode &);
private:
  FitchParsimony &master;
  int column;
  BitSet gapSymbols;
  Symbol gapSymbol;
};


#endif

