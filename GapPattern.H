/****************************************************************
 GapPattern.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_GapPattern_H
#define INCL_GapPattern_H
#include <iostream>
#include "IndelHistory.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;


struct Gap
{
    int pos;
    int len;
    IndelOperation *insertion;
    Gap(int pos=INT_MIN,int len=0,IndelOperation *insertion=NULL)
        : pos(pos), len(len), insertion(insertion) {}
};


class GapPattern
{
public:
    GapPattern();
    int numGaps();
    Gap &operator[](int);
    void appendGap(int pos,int len,IndelOperation * =NULL);
    void printOn(ostream &) const;
    void subsume(GapPattern &other);
    virtual GapPattern *propagateUp(IndelHistory &,int nodeID=-1);
    void sort();
    void coalesce(); // ### OBSOLETE
private:
    Vector<Gap> pattern;
};

ostream &operator<<(ostream &,const GapPattern &);

#endif

