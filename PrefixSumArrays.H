/****************************************************************
 PrefixSumArrays.H
 Bill Majoros - bmajoros@duke.edu
 ****************************************************************/
#ifndef INCL_PrefixSumArrays_H
#define INCL_PrefixSumArrays_H
#include "tigr++/MultiAlignment.H"
#include "tigr++/TigrArray2D.H"


typedef TigrArray2D<int>::RowIn2DArray<int> PrefixSumArrayRef;

class PrefixSumArrays
{
  TigrArray2D<int> matches; // track x alignment pos -> #matches
  TigrArray2D<int> lengths; // track x alignment pos -> cumulative align len

  void init(MultiAlignment &);
public:
  PrefixSumArrays(MultiAlignment &);
  float percentIdentity(int trackNum,int begin,int end);
};


#endif
