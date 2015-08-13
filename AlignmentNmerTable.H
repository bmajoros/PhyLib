/****************************************************************
 AlignmentNmerTable.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_AlignmentNmerTable_H
#define INCL_AlignmentNmerTable_H
#include <iostream>
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/AlphabetMap.H"
#include "BOOM/Vector.H"
#include "BOOM/Array2D.H"
#include "BOOM/Map.H"
#include "BOOM/BitSet.H"
#include "BOOM/Sequence.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                         struct NmerProfile
 ****************************************************************/
struct NmerProfile
{
  int ungappedNmerCode; // only valid for the root
  int gappedNmerCode;
  Vector<int> degenerateMatches; // the non-gapped Nmers that match this one
  Vector<int> degenerateTruncated; // the non-gapped nmers that match this
                                   // one when the last base is ignored
};



/****************************************************************
                      class AlignmentNmerTable
 ****************************************************************/
class AlignmentNmerTable
{
public:
  AlignmentNmerTable(MultSeqAlignment &,AlphabetMap &,int order);
  virtual ~AlignmentNmerTable();
  int getAlignmentLength() const;
  NmerProfile *getNmer(int trackID,int column,int order);
  int getN() const; // =order+1
  AlignmentNmerTable *nextHigherOrder();
  AlignmentNmerTable *getLowerOrder(int order);

protected:
  int alignmentLength; // equals table size +order-1
  int numUngappedNmers, numGappedNmers;
  Array2D<NmerProfile*> table;
  Vector<NmerProfile*> profiles;
  AlphabetMap &alphabetMap;
  Alphabet &gappedAlphabet, &ungappedAlphabet;
  int N; // =order+1
  BOOM::Map<Sequence,NmerProfile*> seqToProfile; // gapped seqs
  const BitSet &gapSymbols;
  MultSeqAlignment &alignment;
  AlignmentNmerTable *lowerOrder;

  void buildTable(MultSeqAlignment &);
  void initDegeneracyLists();
  bool containsGaps(const Sequence &);
};



#endif

