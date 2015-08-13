/****************************************************************
 ResidueCounter.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "ResidueCounter.H"
#include <iostream>
using namespace std;
using namespace BOOM;


void ResidueCounter::increment(char residue)
{
  List<ResidueCount>::iterator cur=countList.begin(),
    end=countList.end();
  for(; cur!=end ; ++cur)
    {
      ResidueCount &count=*cur;
      if(count.residue==residue) {++count.count; return;}
      if(count.residue>residue) 
	{countList.insert(cur,ResidueCount(residue));}
    }
  countList.insert(end,ResidueCount(residue));
}



void ResidueCounter::clear()
{
  countList.clear();
}



void ResidueCounter::add(ResidueCounter &to,ResidueCounter &into)
{
  List<ResidueCount>::iterator 
    aCur=countList.begin(),
    aEnd=countList.end(), 
    bCur=to.countList.begin(),
    bEnd=to.countList.end();
  while(aCur!=aEnd && bCur!=bEnd)
    {
      ResidueCount &a=*aCur, &b=*bCur;
      if(a.residue<b.residue)
	{
	  into.countList.append(a);
	  ++aCur;
	}
      else if(b.residue<a.residue)
	{
	  into.countList.append(b);
	  ++bCur;
	}
      else // a.residue==b.residue
	{
	  into.countList.append(ResidueCount(a.residue,a.count+b.count));
	  ++aCur;
	  ++bCur;
	}
    }
  for(; aCur!=aEnd ; ++aCur) into.countList.append(*aCur);
  for(; bCur!=bEnd ; ++bCur) into.countList.append(*bCur);
}



bool ResidueCounter::getMostPopular(char &residue)
{
  ResidueCount mostPopular(' ');
  mostPopular.count=0;
  bool unique=false;
  List<ResidueCount>::iterator cur=countList.begin(),
    end=countList.end();
  for(; cur!=end ; ++cur)
    {
      ResidueCount &count=*cur;
      if(count.count>mostPopular.count) 
	{
	  mostPopular=count;
	  unique=true;
	}
      else if(count.count==mostPopular.count)
	unique=false;
    }
  if(unique) residue=mostPopular.residue;
  return unique;
}



List<ResidueCount>::iterator ResidueCounter::begin()
{
  return countList.begin();
}



List<ResidueCount>::iterator ResidueCounter::end()
{
  return countList.end();
}




