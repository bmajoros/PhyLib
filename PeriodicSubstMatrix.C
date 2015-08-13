/****************************************************************
 PeriodicSubstMatrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "PeriodicSubstMatrix.H"
using namespace std;
using namespace BOOM;


PeriodicSubstMatrix::PeriodicSubstMatrix(int order,const AlphabetMap &m)
  : matrices(3),
    NthOrdSubstMatrix(order,m)
{
  // ctor
  
  for(int i=0 ; i<3 ; ++i)
    matrices[i]=new NthOrdSubstMatrix(order,m);
}



PeriodicSubstMatrix::~PeriodicSubstMatrix()
{
  for(int i=0 ; i<3 ; ++i)
    delete matrices[i];
}



bool PeriodicSubstMatrix::isPeriodic() const
{
  return true;
}



SubstitutionMatrix *PeriodicSubstMatrix::lookup(const Sequence &context,
						int phase,int begin,int len)
{
  return matrices[phase]->lookup(context,begin,len);
}



SubstitutionMatrix *PeriodicSubstMatrix::lookup(const Sequence &parentContext,
						const Sequence &childContext,
						int phase)
{
  return matrices[phase]->lookup(parentContext,childContext);
}



void PeriodicSubstMatrix::setMatrix(int index,SubstitutionMatrix *M,int phase)
{
  matrices[phase]->setMatrix(index,M);
}



void PeriodicSubstMatrix::setLowerOrderModel(NthOrdSubstMatrix *M,int phase)
{
  matrices[phase]->setLowerOrderModel(M);
}



void PeriodicSubstMatrix::setIthMatrix(int i,NthOrdSubstMatrix *M)
{
  matrices[i]=M;
}



NthOrdSubstMatrix *PeriodicSubstMatrix::getLowerOrderModel(int order,
							   int phase)
{
  return matrices[phase]->getLowerOrderModel(order);
}



SubstitutionMatrix &PeriodicSubstMatrix::getIthMatrix(int i,int phase)
{
  return matrices[phase]->getIthMatrix(i);
}



void PeriodicSubstMatrix::initGenerators()
{
  for(int i=0 ; i<3 ; ++i)
    matrices[i]->initGenerators();
}



Symbol PeriodicSubstMatrix::mutate(Symbol s,
				   const Sequence &context,
				   int contextBegin,
				   int contextLen,
				   int phase)
{
  return matrices[phase]->mutate(s,context,contextBegin,contextLen);
}



void PeriodicSubstMatrix::printOn(ostream &os)
{
  for(int i=0 ; i<3 ; ++i) {
    matrices[i]->printOn(os);
    os<<endl;
  }
}



