/****************************************************************
 IndelHistory.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "IndelHistory.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                      IndelOperation methods
 ****************************************************************/
IndelOperation::IndelOperation(Type type,int length)
    : type(type),
      length(length)
{
    // ctor
}



IndelOperation IndelOperation::operator+(const IndelOperation &other) const
{
    return IndelOperation(type,length+other.length);
}



void IndelOperation::printOn(ostream &os) const
{
    switch(type)
    {
        case INSERTION:
            os<<"i";
            break;
        case DELETION:
            os<<"d";
            break;
        case SUBSTITUTION:
            os<<"s";
            break;
    }
    os<<length;
}



ostream &operator<<(ostream &os,const IndelOperation &op)
{
    op.printOn(os);
    return os;
}



/****************************************************************
                       IndelHistory methods
 ****************************************************************/
IndelHistory::IndelHistory()
{
    // ctor
}



void IndelHistory::append(IndelOperation op)
{
    int L=history.size();
    if(L>0)
    {
        IndelOperation last=history[L-1];
        if(last.getType()==op.getType())
        {
            history[L-1]=last+op;
            return;
        }
    }
    history.push_back(op);
}



void IndelHistory::append(IndelOperation::Type type,int length)
{
    append(IndelOperation(type,length));
}



int IndelHistory::getLength()
{
    return history.size();
}



IndelOperation &IndelHistory::operator[](int i)
{
    return history[i];
}



int IndelHistory::getParentSeqLength()
{
    int L=history.size();
    int effectiveLen=0;
    for(int i=0 ; i<L ; ++i)
    {
        IndelOperation op=history[i];
        if(op.getType()!=IndelOperation::INSERTION)
            effectiveLen+=op.getLength();
    }
    return effectiveLen;
}



int IndelHistory::getChildSeqLength()
{
    int L=history.size();
    int effectiveLen=0;
    for(int i=0 ; i<L ; ++i)
    {
        IndelOperation op=history[i];
        if(op.getType()!=IndelOperation::DELETION)
            effectiveLen+=op.getLength();
    }
    return effectiveLen;
}



void IndelHistory::printOn(ostream &os) const
{
    int n=history.size();
    for(int i=0 ; i<n ; ++i)
        os<<history[i];
}



ostream &operator<<(ostream &os,const IndelHistory &history)
{
    history.printOn(os);
    return os;
}



void IndelHistory::topOff(int sequenceLength)
{
    /*
      This method just extends the history with a final SUBSTITUTION
      op, to take us to the end of the sequence.  Otherwise, the indel
      history would stop at the final insertion or deletion op, due to
      the way indel histories are constructed during evolution.
     */
    
    int pos=0;
    int n=history.size();
    for(int i=0 ; i<n ; ++i)
    {
        IndelOperation op=history[i];
        int opLen=op.getLength();
        switch(op.getType())
        {
            case IndelOperation::INSERTION:    pos+=opLen; break;
            case IndelOperation::DELETION:            break;
            case IndelOperation::SUBSTITUTION: pos+=opLen; break;
        }
    }
    if(pos<sequenceLength)
        append(IndelOperation::SUBSTITUTION,sequenceLength-pos);
}




