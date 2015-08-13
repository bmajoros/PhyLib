/****************************************************************
 DegenerateDnaMatch.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_DegenerateDnaMatch_H
#define INCL_DegenerateDnaMatch_H
#include "BOOM/Sequence.H"
#include "BOOM/BitSet.H"
using namespace std;
using namespace BOOM;


bool degenerateDnaMatch(Sequence &a,Sequence &b,Symbol gapChar,int maxLen=-1); 
bool degenerateDnaMatch(Sequence &a,Sequence &b,const BitSet &gapSymbols,
			int maxLen=-1); 


#endif

