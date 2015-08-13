/****************************************************************
 PeriodicSubstMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_PeriodicSubstMatrix_H
#define INCL_PeriodicSubstMatrix_H
#include <iostream>
#include "NthOrdSubstMatrix.H"
#include "BOOM/Array1D.H"
using namespace std;
using namespace BOOM;


class PeriodicSubstMatrix : public NthOrdSubstMatrix
{
public:
  PeriodicSubstMatrix(int order,const AlphabetMap &);
  virtual ~PeriodicSubstMatrix();
  virtual bool isPeriodic() const;
  virtual SubstitutionMatrix *lookup(const Sequence &context,
				     int phase=-1,
				     int begin=0,int len=-1);
  virtual SubstitutionMatrix *lookup(const Sequence &parentContext,
				     const Sequence &childContext,
				     int phase=-1);
  virtual void setMatrix(int index,SubstitutionMatrix *,int phase=-1);
  virtual void setLowerOrderModel(NthOrdSubstMatrix *,int phase=-1);
  virtual NthOrdSubstMatrix *getLowerOrderModel(int order=-1,int phase=-1);
  virtual SubstitutionMatrix &getIthMatrix(int i,int phase=-1);
  virtual void initGenerators();
  virtual Symbol mutate(Symbol,const Sequence &context,int contextBegin=0,
			int contextLen=-1,int phase=-1);
  virtual void printOn(ostream &);
  void setIthMatrix(int i,NthOrdSubstMatrix*);

private:
  Array1D<NthOrdSubstMatrix*> matrices;
};

#endif

