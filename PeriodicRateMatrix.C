/****************************************************************
 PeriodicRateMatrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "PeriodicRateMatrix.H"
#include "PeriodicSubstMatrix.H"
using namespace std;
using namespace BOOM;


PeriodicRateMatrix::PeriodicRateMatrix(int order,
				       MolecularSequenceType seqType,
				       RateMatrixType matrixType,
				       bool dual) 
  : matrices(3),
    NthOrdRateMatrix(order,seqType,matrixType,dual)
{
  // ctor

  for(int i=0 ; i<3 ; ++i)
    matrices[i]=new NthOrdRateMatrix(order,seqType,matrixType,dual);
}



PeriodicRateMatrix::~PeriodicRateMatrix() 
{
  // dtor

  for(int i=0 ; i<3 ; ++i) delete matrices[i];
}



NthOrdRateMatrix *PeriodicRateMatrix::random(bool dual,
					     double noiseSigma,
					     MolecularSequenceType seqType,
					     int order,
					     RateMatrixType matrixType,
					     BOOM::Array1D<double> 
					           equilibriumFreqs) 
{
  // static method

  PeriodicRateMatrix *M=new PeriodicRateMatrix(order,seqType,matrixType,dual);
  for(int i=0 ; i<3 ; ++i)
    M->matrices[i]=NthOrdRateMatrix::random(dual,noiseSigma,seqType,order,
					    matrixType,equilibriumFreqs);
  return M;
}



/***************************************************************************
  nextHigherOrder() : The recipient of this message becomes useless after
                      this call, so the caller should delete it.  It is not
                      linked into the higher-order matrix.
 */
NthOrdRateMatrix *PeriodicRateMatrix::nextHigherOrder() 
{
  PeriodicRateMatrix *M=new PeriodicRateMatrix(order,seqType,matrixType,dual);
  for(int i=0 ; i<3 ; ++i) {
    M->matrices[i]=matrices[i]->nextHigherOrder();
    matrices[i]=NULL; // to prevent double-deletion
  }
  // delete this;
  return M;
}



RateMatrix &PeriodicRateMatrix::lookup(const Sequence &context,int phase,
				       int begin,int len) 
{
  return matrices[phase]->lookup(context,begin,len);
}



RateMatrix &PeriodicRateMatrix::lookup(const Sequence &parentContext,
				       const Sequence &childContext,
				       int phase) 
{
  return matrices[phase]->lookup(parentContext,childContext,phase);
}



RateMatrix &PeriodicRateMatrix::getIthMatrix(int i,int phase) 
{
  return matrices[phase]->getIthMatrix(i);
}



NthOrdSubstMatrix *PeriodicRateMatrix::instantiate(double branchLength) 
{
  PeriodicSubstMatrix *M=new PeriodicSubstMatrix(order,alphabetMap);
  for(int i=0 ; i<3 ; ++i)
    M->setIthMatrix(i,matrices[i]->instantiate(branchLength));
  return M;
}



NthOrdRateMatrix *PeriodicRateMatrix::getLowerOrderModel(int order,
							 int phase) 
{
  return matrices[phase]->getLowerOrderModel(order);
}



NthOrdRateMatrix *PeriodicRateMatrix::ancestralContextsOnly() 
{
  PeriodicRateMatrix *M=new PeriodicRateMatrix(order,seqType,matrixType,dual);
  for(int i=0 ; i<3 ; ++i)
    M->matrices[i]=matrices[i]->ancestralContextsOnly();
  return M;
}



NthOrdRateMatrix *PeriodicRateMatrix::clone() 
{
  PeriodicRateMatrix *M=new PeriodicRateMatrix(order,seqType,matrixType,dual);
  for(int i=0 ; i<3 ; ++i)
    M->matrices[i]=matrices[i]->clone();
  return M;
}



bool PeriodicRateMatrix::isPeriodic() 
{
  return true;
}



void PeriodicRateMatrix::save(const BOOM::String &filename) 
{
  ofstream os(filename.c_str());
  if(!os.good()) throw filename+" : can't write to file";
  save(os);
}



void PeriodicRateMatrix::save(ostream &os) 
{
  os<<1; // 1=periodic (yes)
  for(int i=0 ; i<3 ; ++i)
    matrices[i]->save(os);
}



void PeriodicRateMatrix::load_noIndicator(istream &is) 
{
  for(int i=0 ; i<3 ; ++i)
    matrices[i]=NthOrdRateMatrix::load(is);
  NthOrdRateMatrix *M=matrices[0];
  order=M->getOrder();
  seqType=M->getSeqType();
  matrixType=M->getMatrixType();
  dual=M->isDual();
}



