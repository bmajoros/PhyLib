/****************************************************************
 DegenerateDnaMatch.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "DegenerateDnaMatch.H"
using namespace std;
using namespace BOOM;


bool degenerateDnaMatch(Sequence &a,Sequence &b,Symbol gapChar,int maxLen) 
{
  int n=(maxLen<0 ? a.getLength() : maxLen);
  for(int i=0 ; i<n ; ++i) {
    Symbol sa=a[i], sb=b[i];
    if(sa!=gapChar && sb!=gapChar && sa!=sb)
      return false;
  }
  return true;
}



bool degenerateDnaMatch(Sequence &a,Sequence &b,const BitSet &gapSymbols,
			int maxLen) 
{
  int n=(maxLen<0 ? a.getLength() : maxLen);
  for(int i=0 ; i<n ; ++i) {
    Symbol sa=a[i], sb=b[i];
    if(!gapSymbols.isMember(sa) && !gapSymbols.isMember(sb) && sa!=sb)
      return false;
  }
  return true;
}

