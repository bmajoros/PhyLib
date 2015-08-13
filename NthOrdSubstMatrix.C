/****************************************************************
 NthOrdSubstMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "NthOrdSubstMatrix.H"
using namespace BOOM;


NthOrdSubstMatrix::NthOrdSubstMatrix(const NthOrdSubstMatrix &other)
  : alphabetMap(other.alphabetMap), order(other.order),
    matrices(other.matrices.size()), lowerOrderModel(NULL)
{
  // copy ctor

  matrices.setAllTo(NULL);
  int n=matrices.size();
  for(int i=0 ; i<n ; ++i)
    if(other.matrices[i])
      matrices[i]=new SubstitutionMatrix(*other.matrices[i]);
  if(other.lowerOrderModel)
    lowerOrderModel=new NthOrdSubstMatrix(*other.lowerOrderModel);
}



NthOrdSubstMatrix::NthOrdSubstMatrix(int order,
				     const AlphabetMap &alphabetMap)
  : lowerOrderModel(NULL),
    alphabetMap(alphabetMap),
    order(order)
{
  // ctor

  float radix=alphabetMap.getRangeSize();
  int numMatrices=pow(radix,float(order));
  matrices.resize(numMatrices);
  matrices.setAllTo(NULL);
}



NthOrdSubstMatrix::~NthOrdSubstMatrix()
{
  // ctor

  int n=matrices.size();
  for(int i=0 ; i<n ; ++i)
    delete matrices[i];

  delete lowerOrderModel;
}



SubstitutionMatrix *NthOrdSubstMatrix::lookup(unsigned contextCode,
					      int contextLength)
{
  if(contextLength==order)
    return matrices[contextCode];
  else 
    return lowerOrderModel->lookup(contextCode,contextLength);
}



SubstitutionMatrix *NthOrdSubstMatrix::lookup(const Sequence &context,
					      int phase,int begin,int len)
{
  int L=len;
  if(L<0) L=context.getLength();
  if(L==order) {
    unsigned index=context.asInt(alphabetMap,begin,len);
    return matrices[index];
  }
  else if(L<order)
    return lowerOrderModel->lookup(context,begin,len);
  else {
      cout<<"begin="<<begin<<", len="<<len<<", L="<<L<<", context.len="<<context.getLength()<<", order="<<order<<endl;
      throw "Bad context in NthOrdSubstMatrix::lookup()";
  }
}



SubstitutionMatrix *NthOrdSubstMatrix::lookup(const Sequence &parentContext,
                                              const Sequence &childContext,
					      int phase) {
    Sequence context=parentContext;
    context.append(childContext);
    return lookup(context);
}



void NthOrdSubstMatrix::setMatrix(int index,SubstitutionMatrix *M,int phase)
{
  matrices[index]=M;
}



void NthOrdSubstMatrix::setLowerOrderModel(NthOrdSubstMatrix *M,int phase)
{
  lowerOrderModel=M;
}



int NthOrdSubstMatrix::getNumMatrices() const {
    return matrices.size();
}



SubstitutionMatrix &NthOrdSubstMatrix::getIthMatrix(int i,int phase) {
  return *matrices[i];
}



int NthOrdSubstMatrix::getOrder() const {
    return order;
}



NthOrdSubstMatrix *NthOrdSubstMatrix::getLowerOrderModel(int order,
							 int phase) {
    if(order<0 || order==this->order-1) return lowerOrderModel;
    return lowerOrderModel->getLowerOrderModel(order);
}



void NthOrdSubstMatrix::initGenerators() {
    int n=matrices.size();
    for(int i=0 ; i<n ; ++i) {
        SubstitutionMatrix *M=matrices[i];
        if(M) M->initGenerators();
    }
    if(lowerOrderModel) lowerOrderModel->initGenerators();
}



Symbol NthOrdSubstMatrix::mutate(Symbol s,const Sequence &context,
                                 int contextBegin,int contextLen,int phase)
{
    SubstitutionMatrix *M=lookup(context,contextBegin,contextLen);
    return M->mutate(s);
}



void NthOrdSubstMatrix::printOn(ostream &os) {
    int n=matrices.size();
    for(int i=0 ; i<n ; ++i) {
        os<<*matrices[i]<<endl;
    }
}



ostream &operator<<(ostream &os,NthOrdSubstMatrix &M) {
    M.printOn(os);
    return os;
}



bool NthOrdSubstMatrix::isPeriodic() const
{
  return false;
}



void NthOrdSubstMatrix::convertToLogs() 
{
  int n=matrices.size();
  for(int i=0 ; i<n ; ++i) {
    SubstitutionMatrix *M=matrices[i];
    if(M) M->convertToLogs();
  }
}



NthOrdSubstMatrix *NthOrdSubstMatrix::clone() const
{
  return new NthOrdSubstMatrix(*this);
}
