/****************************************************************
 ResidueCounter.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ResidueCounter_H
#define INCL_ResidueCounter_H
#include "BOOM/List.H"
using namespace BOOM;

/****************************************************************
                       struct ResidueCount
 ****************************************************************/
struct ResidueCount
{
  char residue;
  int count;
  ResidueCount(char r,int c=1) : residue(r), count(c) {}
};


/****************************************************************
                       class ResidueCounter
 ****************************************************************/
class ResidueCounter
{
  BOOM::List<ResidueCount> countList;
public:
  void increment(char residue);
  void clear();
  void add(ResidueCounter &to,ResidueCounter &into);
  bool getMostPopular(char &residue);
  BOOM::List<ResidueCount>::iterator begin();
  BOOM::List<ResidueCount>::iterator end();
};


#endif
