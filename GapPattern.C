/****************************************************************
 GapPattern.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "GapPattern.H"
#include "BOOM/VectorSorter.H"
using namespace std;
using namespace BOOM;


GapPattern::GapPattern()
{
    // ctor
}



int GapPattern::numGaps()
{
    return pattern.size();
}



Gap &GapPattern::operator[](int i)
{
    return pattern[i];
}



void GapPattern::appendGap(int pos,int len,IndelOperation *insertion)
{
    pattern.push_back(Gap(pos,len,insertion));
}



void GapPattern::printOn(ostream &os) const
{
    int n=pattern.size();
    for(int i=0 ; i<n ; ++i)
    {
        Gap gap=pattern[i];
        os<<"("<<gap.pos<<","<<gap.len
          <<") ";
    }
}



ostream &operator<<(ostream &os,const GapPattern &pattern)
{
    pattern.printOn(os);
    return os;
}



void GapPattern::subsume(GapPattern &other)
{
    int n1=pattern.size(), n2=other.pattern.size();
    int i=0;
    for(int j=0 ; j<n2 ; ++j)
    {
        Gap otherGap=other[j];
        while(i<n1 && pattern[i].pos<=otherGap.pos) ++i;
        pattern.insertByIndex(otherGap,i);
        ++n1;
    }
}



GapPattern *GapPattern::propagateUp(IndelHistory &h,int nodeID)
{
    IndelHistory H=h;
    appendGap(H.getChildSeqLength(),0); // to force full examination of H
    GapPattern &newPattern=*new GapPattern;
    int nH=H.getLength(), nG=pattern.size();
    int parentPos=0, childPos=0;
    for(int iH=0, iG=0 ; iG<nG ; ++iG)
    {
        Gap oldGap=pattern[iG];
        int newGapPos;
        for( ; iH<nH ; ++iH)
        {
            IndelOperation &op=H[iH];
            int nextParentPos, nextChildPos, opLen=op.getLength();
            switch(op.getType())
            {
                case IndelOperation::SUBSTITUTION:
                    nextParentPos=parentPos+opLen;
                    nextChildPos=childPos+opLen;
                    if(childPos<=oldGap.pos && nextChildPos>oldGap.pos) {
                        newGapPos=parentPos+oldGap.pos-childPos;
                        goto ENDLOOP;
                    }
                    break;
                case IndelOperation::INSERTION:
                    nextParentPos=parentPos;
                    nextChildPos=childPos+opLen;
                    if(childPos<=oldGap.pos && nextChildPos>oldGap.pos)
                    { // "insertion within an insertion"
                        int leftPart=oldGap.pos-childPos;
                        if(leftPart>0) {
                            newPattern.appendGap(parentPos,leftPart,&h[iH]);
                            op.changeLength(opLen-leftPart);
                            childPos+=leftPart;
                        }
                        newGapPos=parentPos;
                        goto ENDLOOP;
                    }
                    else newPattern.appendGap(parentPos,opLen,&h[iH]);
                    break;
                case IndelOperation::DELETION:
                    nextParentPos=parentPos+opLen;
                    nextChildPos=childPos;
                    break;
            }
            parentPos=nextParentPos;
            childPos=nextChildPos;
        }
      ENDLOOP:
        if(oldGap.len>0)
            newPattern.appendGap(newGapPos,oldGap.len,oldGap.insertion);
    }
    return &newPattern;
}



class GapComparator : public Comparator<Gap>
{
public:
    virtual bool equal(Gap &a,Gap &b)   {return a.pos==b.pos;}
    virtual bool greater(Gap &a,Gap &b) {return a.pos>b.pos;}
    virtual bool less(Gap &a,Gap &b)    {return a.pos<b.pos;}
};



void GapPattern::sort()
{
    GapComparator cmp;
    VectorSorter<Gap> sorter(pattern,cmp);
    sorter.sortAscendInPlace();
}



void GapPattern::coalesce()
{
    throw "GapPattern::coalesce() -- shouldn't be calling this!";
    
    int n=pattern.size()-1;
    for(int i=0 ; i<n ; ++i)
    {
        Gap &a=pattern[i], &b=pattern[i+1];
        if(a.pos==b.pos)
        {
            a.len+=b.len;
            pattern.cut(i+1,1);
            --n;
            --i;
        }
    }
}



