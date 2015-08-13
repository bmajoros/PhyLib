/****************************************************************
 NthOrdSubstMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_NthOrdSubstMatrix_H
#define INCL_NthOrdSubstMatrix_H
#include <iostream>
#include "SubstitutionMatrix.H"
#include "BOOM/Array1D.H"
#include "BOOM/Sequence.H"
#include "BOOM/Symbol.H"
using namespace std;
using namespace BOOM;

class NthOrdSubstMatrix
{
public:
  // Constructors/destructors
  NthOrdSubstMatrix(int order,const AlphabetMap &);
  NthOrdSubstMatrix(const NthOrdSubstMatrix &);
  virtual ~NthOrdSubstMatrix();

  // Instance methods
  int getNumMatrices() const;
  int getOrder() const;
  NthOrdSubstMatrix *clone() const;

  // Virtual methods
  virtual void convertToLogs();
  virtual bool isPeriodic() const;
  virtual SubstitutionMatrix *lookup(unsigned contextCode,int contextLength);
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
protected:
  Array1D<SubstitutionMatrix*> matrices;
  NthOrdSubstMatrix *lowerOrderModel;
  int order;
  const AlphabetMap &alphabetMap;
};

ostream &operator<<(ostream &,NthOrdSubstMatrix &);

#endif
