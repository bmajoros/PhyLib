/****************************************************************
 PrefixSumArrays.C
 Bill Majoros - bmajoros@duke.edu
 ****************************************************************/
#include "PrefixSumArrays.H"
#include <iostream>
using namespace std;


PrefixSumArrays::PrefixSumArrays(MultiAlignment &alignment)
  : matches(alignment.getNumTracks()-1,alignment.getLength()),
    lengths(alignment.getNumTracks()-1,alignment.getLength())
{
  // ctor

  init(alignment);
}



float PrefixSumArrays::percentIdentity(int trackNum,int begin,int end)
{
  int trackNumMinus1=trackNum-1;
  PrefixSumArrayRef matchesArray=matches[trackNumMinus1];
  PrefixSumArrayRef lengthArray=lengths[trackNumMinus1];
  int numMatches=matchesArray[end]-matchesArray[begin];
  int length=lengthArray[end]-lengthArray[begin];
  float percent=numMatches/float(length);
  return percent;
}



void PrefixSumArrays::init(MultiAlignment &alignment)
{
  AlignmentTrack &targetTrack=alignment.getIthTrack(0);
  int numTracks=alignment.getNumTracks();
  int alignmentLength=alignment.getLength();
  for(int trackNum=1 ; trackNum<numTracks ; ++trackNum)
    {
      int M=0, L=0;
      AlignmentTrack &track=alignment.getIthTrack(trackNum);
      int trackLen=track.getLength();
      PrefixSumArrayRef matchArray=matches[trackNum-1];
      PrefixSumArrayRef lengthArray=lengths[trackNum-1];
      for(int pos=0 ; pos<alignmentLength ; ++pos)
	{
	  if(pos<trackLen) 
	    {
	      char a=targetTrack[pos];
	      char b=track[pos];
	      if(a==b)
		if(a!='-') {++M; ++L;}
		else {} // both dashes -- do nothing
	      else ++L; // mismatch or indel
	    }
	  matchArray[pos]=M;
	  lengthArray[pos]=L;
	}
    }
}


