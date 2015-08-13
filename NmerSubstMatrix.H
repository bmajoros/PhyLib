/****************************************************************
 NmerSubstMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_NmerSubstMatrix_H
#define INCL_NmerSubstMatrix_H
#include "NthOrdSubstMatrix.H"
using namespace std;
using namespace BOOM;


class NmerSubstMatrix : public NthOrdSubstMatrix
{
public:
  NmerSubstMatrix(int order,const AlphabetMap &);
  double lookupNmers(const Sequence &parentNmer,const Sequence &childNmer);
  double lookupNmers(int parentNmerCode,int childNmerCode);
  virtual void printOn(ostream &);
  virtual void convertToLogs();
  GSL::Matrix &peek();
  void averageLowerOrderModel(const Array1D<double> &);
protected:
  GSL::Matrix M;
  int numNmers;
};

#endif

