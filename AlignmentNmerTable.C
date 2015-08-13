/****************************************************************
 AlignmentNmerTable.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "math.h"
#include "AlignmentNmerTable.H"
#include "DegenerateDnaMatch.H"
using namespace std;
using namespace BOOM;


AlignmentNmerTable::AlignmentNmerTable(MultSeqAlignment &alignment,
				       AlphabetMap &alphabetMap,
				       int order)
  : alphabetMap(alphabetMap),
    N(order+1),
    alignmentLength(alignment.getLength()),
    numUngappedNmers(pow(alphabetMap.getRangeSize(),order+1)),
    numGappedNmers(pow(alphabetMap.getDomainSize(),order+1)),
    table(alignment.getNumTracks(),alignment.getLength()-order),
    gappedAlphabet(*alphabetMap.getDomain()),
    ungappedAlphabet(*alphabetMap.getRange()),
    gapSymbols(alignment.getGapSymbols()),
    alignment(alignment),
    lowerOrder(NULL)
{
  // ctor

  buildTable(alignment);
}



AlignmentNmerTable::~AlignmentNmerTable()
{
  // dtor

  Vector<NmerProfile*>::iterator cur=profiles.begin(), end=profiles.end();
  for(; cur!=end ; ++cur) {
    NmerProfile *victim=*cur;
    delete victim;
  }
}



void AlignmentNmerTable::buildTable(MultSeqAlignment &alignment)
{
  int tableLen=alignment.getLength()-N+1;
  Sequence gappedNmer, ungappedNmer;
  int numTracks=alignment.getNumTracks();
  for(int trackID=0 ; trackID<numTracks ; ++trackID) {
    AlignmentSeq &track=alignment.getIthTrack(trackID);
    for(int col=0 ; col<tableLen ; ++col) {
      const Sequence &trackSeq=track.getSeq();
      if(trackSeq.getLength()<col+N) continue;
      trackSeq.getSubsequence(col,N,gappedNmer);
      NmerProfile *profile;
      if(seqToProfile.isDefined(gappedNmer)) {
	profile=seqToProfile[gappedNmer];
      }
      else {
	seqToProfile[gappedNmer]=profile=new NmerProfile;
	profiles.push_back(profile);
	profile->gappedNmerCode=gappedNmer.asInt(gappedAlphabet);
	gappedNmer.translate(alphabetMap,ungappedNmer);
	profile->ungappedNmerCode=ungappedNmer.asInt(ungappedAlphabet);
      }
      table[trackID][col]=profile;
    }
  }
  initDegeneracyLists();
  seqToProfile.clear();
}



int AlignmentNmerTable::getAlignmentLength() const
{
  return alignmentLength;
}



int AlignmentNmerTable::getN() const
{
  return N;
}



NmerProfile *AlignmentNmerTable::getNmer(int trackID,int column,int order)
{
  if(order<N-1) return lowerOrder->getNmer(trackID,column,order);
  return table[trackID][column];
}



void AlignmentNmerTable::initDegeneracyLists()
{
  Sequence nmer1, nmer2, ungappedNmer; //, truncated;
  int n=profiles.size();
  Vector<NmerProfile*>::iterator cur=profiles.begin(), end=profiles.end();
  for(; cur!=end ; ++cur) {
    NmerProfile *profile=*cur;
    nmer1.fromInt(profile->gappedNmerCode,N,gappedAlphabet);
    for(int nmerCode=0 ; nmerCode<numGappedNmers ; ++nmerCode) {
      nmer2.fromInt(nmerCode,N,gappedAlphabet);
      if(containsGaps(nmer2)) continue;
      if(degenerateDnaMatch(nmer1,nmer2,gapSymbols)) {
	nmer2.translate(alphabetMap,ungappedNmer);
	int ungappedNmerCode=ungappedNmer.asInt(ungappedAlphabet);
	profile->degenerateMatches.push_back(ungappedNmerCode);
      }
      if(degenerateDnaMatch(nmer1,nmer2,gapSymbols,N-1)) {
	nmer2.translate(alphabetMap,ungappedNmer);
	int ungappedNmerCode=ungappedNmer.asInt(ungappedAlphabet);
	profile->degenerateTruncated.push_back(ungappedNmerCode);
      }
    /*
      if(N>1) {
	nmer2.getSubsequence(0,N-1,truncated);
	int code=truncated.asInt(ungappedAlphabet);
	profile->degenerateTruncated.push_back(code);
      }
    */
    }
  }
}



bool AlignmentNmerTable::containsGaps(const Sequence &seq)
{
  int L=seq.getLength();
  for(int i=0 ; i<L ; ++i)
    if(gapSymbols.isMember(seq[i]))
      return true;
  return false;
}



AlignmentNmerTable *AlignmentNmerTable::getLowerOrder(int order)
{
  return lowerOrder;
}



AlignmentNmerTable *AlignmentNmerTable::nextHigherOrder()
{
  AlignmentNmerTable *t=new AlignmentNmerTable(alignment,alphabetMap,N);
  t->lowerOrder=this;
  return t;
}



