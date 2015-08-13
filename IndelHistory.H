/****************************************************************
 IndelHistory.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_IndelHistory_H
#define INCL_IndelHistory_H
#include <iostream>
#include "BOOM/Vector.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                        class IndelOperation
 ****************************************************************/
class IndelOperation
{
public:
    enum Type { INSERTION, DELETION, SUBSTITUTION };
    IndelOperation(Type,int length=1);
    inline Type getType();
    inline int getLength();
    inline void extendBy1();
    IndelOperation operator+(const IndelOperation &) const;
    void printOn(ostream &) const;
    void changeLength(int l) {length=l;}
private:
    Type type;
    int length;
};

ostream &operator<<(ostream &,const IndelOperation &);



/****************************************************************
                        class IndelHistory
 ****************************************************************/
class IndelHistory
{
public:
    IndelHistory();
    void append(IndelOperation);
    void append(IndelOperation::Type,int length=1);
    int getLength();
    int getParentSeqLength();
    int getChildSeqLength();
    IndelOperation &operator[](int);
    void printOn(ostream &) const;
    void topOff(int sequenceLength);
private:
    Vector<IndelOperation> history;
};

ostream &operator<<(ostream &,const IndelHistory &);



/****************************************************************
                         inline methods
 ****************************************************************/
IndelOperation::Type IndelOperation::getType()
{
    return type;
}



int IndelOperation::getLength()
{
    return length;
}



void IndelOperation::extendBy1()
{
    ++length;
}



#endif

