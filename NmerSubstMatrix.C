/****************************************************************
 NmerSubstMatrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "NmerSubstMatrix.H"
#include "BOOM/DnaAlphabet.H"
using namespace std;
using namespace BOOM;


NmerSubstMatrix::NmerSubstMatrix(int order,const AlphabetMap &alphabetMap)
  : NthOrdSubstMatrix(order,alphabetMap)
{
  // ctor

  int n=order+1;
  numNmers=pow((float)alphabetMap.getRangeSize(),(float)n);
  M.resize(numNmers,numNmers);

  if(order>0) {
    delete lowerOrderModel;
    lowerOrderModel=new NmerSubstMatrix(order-1,alphabetMap);
  }
}



double NmerSubstMatrix::lookupNmers(int parentCode,int childCode)
{
  return M(parentCode,childCode);
}



double NmerSubstMatrix::lookupNmers(const Sequence &parentNmer,
				    const Sequence &childNmer)
{
  int len=parentNmer.getLength();
  if(len<=order) {
    NmerSubstMatrix *M=static_cast<NmerSubstMatrix*>(lowerOrderModel);
    return M->lookupNmers(parentNmer,childNmer);
  }
  int parentCode=parentNmer.asInt(alphabetMap);
  int childCode=childNmer.asInt(alphabetMap);

  return M(parentCode,childCode);
}



void NmerSubstMatrix::printOn(ostream &os)
{
  for(int i=0 ; i<numNmers ; ++i) {
    double sum=0;
    for(int j=0 ; j<numNmers ; ++j) {
      double x=exp(M(i,j));
      os<<x<<"\t";
      sum+=x;
    }
    os<<" SUM["<<i<<"]="<<sum<<endl;
  }
  os<<endl;
      
  os<<M;
}



GSL::Matrix &NmerSubstMatrix::peek()
{
  return M;
}



void NmerSubstMatrix::averageLowerOrderModel(const Array1D<double> &eq)
{
  if(!lowerOrderModel) return;
  Alphabet &myAlphabet=*alphabetMap.getRange();
  int numAlpha=myAlphabet.size();
  NmerSubstMatrix *lowerModel=static_cast<NmerSubstMatrix*>(lowerOrderModel);
  GSL::Matrix &lowerM=lowerModel->M;
  int n=order+1, lowerN=order;
  int numLowerNmers=lowerModel->numNmers;
  Sequence lowFrom, lowTo;
  for(int i=0 ; i<numLowerNmers ; ++i) {
    lowFrom.fromInt(i,lowerN,myAlphabet);
    for(int j=0 ; j<numLowerNmers ; ++j) {
      lowTo.fromInt(j,lowerN,myAlphabet);
      double &lowerEntry=lowerM(i,j);
      lowerEntry=0;
      for(Symbol x=0 ; x<numAlpha ; ++x) {
	Sequence hiFrom=lowFrom+x;
	int hiFromIndex=hiFrom.asInt(alphabetMap);
	for(Symbol y=0 ; y<numAlpha ; ++y) {
	  Sequence hiTo=lowTo+y;
	  int hiToIndex=hiTo.asInt(alphabetMap);
	  double hiEntry=M(hiFromIndex,hiToIndex);
	  lowerEntry+=hiEntry*eq[x];
	}
      }
    }
  }

  lowerModel->averageLowerOrderModel(eq);
}



void NmerSubstMatrix::convertToLogs()
{
  for(int i=0 ; i<numNmers ; ++i)
    for(int j=0 ; j<numNmers ; ++j)
      M(i,j)=log(M(i,j));
}


