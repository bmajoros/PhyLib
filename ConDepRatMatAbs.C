/****************************************************************
 ConDepRatMatAbs.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "ConDepRatMatAbs.H"
#include "BOOM/GSL/LabeledMatrixLoader.H"
#include "BOOM/GSL/Vector.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/AminoAlphabet.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/Random.H"
#include "BOOM/Sequence.H"
#include "BOOM/NthOrderStringIterator.H"
#include "BOOM/GSL/GaussianDistribution.H"
#include "JukesCantor.H"
#include "FEL.H"
#include "Kimura2Param.H"
#include "HKY.H"
#include "REV.H"
#include <fstream>
using namespace BOOM;


BOOM::Array1D<double> ConDepRatMatAbs::omitted;

/*
  NOTE: In dual models (those for which we model both parental and child
  contexts), since the parent & child context lengths are assumed equal,
  there can only be even orders.  However, we keep the odd orders anyway,
  because they can be used to generate successively higher-order models 
  via Gaussian perturbation (for simulation runs).
 */

ConDepRatMatAbs::ConDepRatMatAbs(int order,MolecularSequenceType seqType,
				   RateMatrixType matrixType,bool dual)
  : order(order), 
    seqType(seqType), 
    matrixType(matrixType),
    lowerOrderModel(NULL),
    alphabet(seqType==DNA ? (Alphabet&) DnaAlphabet::global : 
             (Alphabet&) AminoAlphabet::global),
    dual(dual)
{
  // ctor

  if(seqType==DNA)
    alphabetMap=DropGapMapping(alphabet,PureDnaAlphabet::global());
  else
    alphabetMap=AlphabetIdentityMap(alphabet);

  int trueOrder=order;
  if(dual && order%2>0) trueOrder=order/2+1; // ODD ORDER IN DUAL MODEL
  int numMatrices=pow(float(alphabetMap.getRangeSize()),float(trueOrder));
  matrices.resize(numMatrices);
  matrices.setAllTo(NULL);

  if(order>0) 
      lowerOrderModel=new ConDepRatMatAbs(order-1,seqType,matrixType,dual);
}



ConDepRatMatAbs::~ConDepRatMatAbs()
{
  int n=matrices.size();
  for(int i=0 ; i<n ; ++i)
    delete matrices[i];

  delete lowerOrderModel;
}



void ConDepRatMatAbs::randomize(double noiseSigma,
                                 BOOM::Array1D<double> equilibriumFreqs) {

    // PRECONDITION: the model is NOT a dual model
    
    GSL::GaussianDistribution noiseModel(0,noiseSigma);
    if(!lowerOrderModel)
        matrices[0]=RateMatrix::random(seqType,matrixType,equilibriumFreqs);
    else {
        lowerOrderModel->randomize(noiseSigma,equilibriumFreqs);
        Array1D<RateMatrix*> &lowerOrderMatrices=lowerOrderModel->matrices;
        Sequence nmer;
        int numMatrices=matrices.size();
        for(int i=0 ; i<numMatrices ; ++i) {
            nmer.fromInt(i,order,alphabetMap);
            int lowerOrderIndex=nmer.asInt(alphabetMap,1,order-1);
            RateMatrix *M=lowerOrderMatrices[lowerOrderIndex];
            RateMatrix *newM=M->clone();
            newM->addNoise(noiseModel);
            matrices[i]=newM;
        }
    }
}



void ConDepRatMatAbs::randomizeDual(double noiseSigma,
                                     BOOM::Array1D<double> equilibriumFreqs) {
    Sequence nmer, parentNmer, childNmer;
    GSL::GaussianDistribution noiseModel(0,noiseSigma);

    if(!lowerOrderModel) // ZEROTH ORDER
        matrices[0]=RateMatrix::random(seqType,matrixType,equilibriumFreqs);

    else if(order%2>0) { // ODD ORDER => PARENT CONTEXT ONLY
        lowerOrderModel->randomizeDual(noiseSigma,equilibriumFreqs);
        ConDepRatMatAbs *lowerOddOrderModel=lowerOrderModel;
        if(lowerOrderModel->lowerOrderModel)
            lowerOddOrderModel=lowerOrderModel->lowerOrderModel;
        Array1D<RateMatrix*> &lowerOrderMatrices=lowerOddOrderModel->matrices;
        int trueOrder=order/2+1;
        int numNmers=pow((float)alphabetMap.getRangeSize(),(float)trueOrder);
        for(int nmerCode=0 ; nmerCode<numNmers ; ++nmerCode) {
            parentNmer.fromInt(nmerCode,trueOrder,alphabetMap);
            int suffixCode=parentNmer.asInt(alphabetMap,1,trueOrder-1);
            RateMatrix *baseModel=lowerOrderMatrices[suffixCode];
            RateMatrix *newModel=baseModel->clone();
            newModel->addNoise(noiseModel);
            matrices[nmerCode]=newModel;
        }
    }

    else { // EVEN ORDER (GREATER THAN ZERO) => PARENT & CHILD CONTEXTS
        int halfOrder=order/2;
        lowerOrderModel->randomizeDual(noiseSigma,equilibriumFreqs);
        Array1D<RateMatrix*> &lowerOrderMatrices=lowerOrderModel->matrices;
        int numMatrices=matrices.size();
        for(int i=0 ; i<numMatrices ; ++i) {
            nmer.fromInt(i,order,alphabetMap);
            int parentNmerCode=nmer.asInt(alphabetMap,0,halfOrder);
            RateMatrix *baseModel=lowerOrderMatrices[parentNmerCode];
            RateMatrix *newM=baseModel->clone();
            newM->addNoise(noiseModel);
            matrices[i]=newM;
        }
    }
}



ConDepRatMatAbs *ConDepRatMatAbs::random(bool dual,
                                           double noiseSigma,
                                           MolecularSequenceType seqType,
                                           int order,
                                           RateMatrixType matrixType,
                                           BOOM::Array1D<double> 
                                           equilibriumFreqs=omitted) {
    ConDepRatMatAbs *Q=new ConDepRatMatAbs(order,seqType,matrixType,dual);
    if(dual)
        Q->randomizeDual(noiseSigma,equilibriumFreqs);
    else
        Q->randomize(noiseSigma,equilibriumFreqs);
    return Q;
}



void ConDepRatMatAbs::save(const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  save(os);
}



void ConDepRatMatAbs::save(ostream &os)
{
    os<<(dual?1:0)<<endl;
    os<<seqType<<endl;
    os<<matrixType<<endl;
    os<<order<<endl;
    unsigned n=matrices.size();
    os<<n<<endl;
    RateMatrix *Q0=order>0 ? getLowerOrderModel(0)->matrices[0] : NULL;
    for(unsigned i=0 ; i<n ; ++i) {
        RateMatrix *Q=matrices[i];
	if(!Q) Q=Q0; // odd-order models usually are not allocated...
        Q->save(os);
    }
    
    if(lowerOrderModel)
        lowerOrderModel->save(os);
}



ConDepRatMatAbs *ConDepRatMatAbs::load(const BOOM::String &filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw BOOM::String("Error opening file: ")+filename;
  return load(is);
}



ConDepRatMatAbs *ConDepRatMatAbs::load(istream &is)
{
  int numMatrices;
  MolecularSequenceType seqType;
  RateMatrixType matrixType;
  int order, dualInt;
  is>>dualInt;
  bool dual=dualInt;
  is>>seqType>>matrixType>>order>>numMatrices;
  ConDepRatMatAbs *Q=new ConDepRatMatAbs(order,seqType,matrixType,dual);
  for(int i=0 ; i<numMatrices ; ++i)
    Q->matrices[i]=RateMatrix::load(is);

  if(order>0)
    Q->lowerOrderModel=load(is);

  return Q;
}



RateMatrixType ConDepRatMatAbs::getType()
{
  return matrixType;
}



int ConDepRatMatAbs::getOrder() const
{
  return order;
}



RateMatrix &ConDepRatMatAbs::lookup(const Sequence &context,int begin,
				     int len)
{
    int L=len;
    if(L<0) L=context.getLength();
    if(L==order) {
        int index=context.asInt(alphabetMap,begin,len);
        return *matrices[index];
    }
    else if(L<order)
        return lowerOrderModel->lookup(context,begin,len);
    else
        throw "Invalid context in ConDepRatMatAbs::lookup()";
}



int ConDepRatMatAbs::lookupIndex(const Sequence &parentContext,
                                  const Sequence &childContext) {
    // PRECONDITION: this is a DUAL model
    
    int L=parentContext.getLength()+childContext.getLength();
    if(L==order) {
        Sequence context=parentContext;
        context.append(childContext);
        return context.asInt(alphabetMap);
    }
    else {
        cout<<"parentContext:";
        parentContext.printOn(cout,alphabet);
        cout<<"\nchildContext:";
        childContext.printOn(cout,alphabet);
        cout<<endl;
        throw "Invalid context in ConDepRatMatAbs::lookupIndex()";
    }
}



RateMatrix &ConDepRatMatAbs::lookup(const Sequence &parentContext,
                                     const Sequence &childContext) {
    Sequence context=parentContext;
    context.append(childContext);
    return lookup(context);
}



NthOrdSubstMatrix *ConDepRatMatAbs::instantiate(double branchLength)
{
  NthOrdSubstMatrix *Pt=new NthOrdSubstMatrix(order,alphabetMap);
  instantiate(branchLength,*Pt);
  return Pt;
}



void ConDepRatMatAbs::instantiate(double branchLength,NthOrdSubstMatrix &Pt)
{
    unsigned numMatrices=matrices.size();
    for(unsigned i=0 ; i<numMatrices ; ++i) {
        RateMatrix *M=matrices[i];
        if(!M) continue;
        SubstitutionMatrix *matrix=M->instantiate(branchLength);
        Pt.setMatrix(i,matrix);
    }
    if(lowerOrderModel) {
        Pt.setLowerOrderModel(lowerOrderModel->instantiate(branchLength));
    }
}



int ConDepRatMatAbs::getNumRows() const
{
  RateMatrix *Q=matrices[0];
  if(!Q) throw "NULL matrix in ConDepRatMatAbs::getNumRows()";
  return Q->getNumRows();
}



const Alphabet &ConDepRatMatAbs::getAlphabet() const
{
  return alphabet;
}



int ConDepRatMatAbs::numParameters() const
{
  RateMatrix *Q=matrices[0];
  if(!Q) throw "NULL matrix in ConDepRatMatAbs::numParameters()";
  return Q->numParameters();
}



RateMatrixType ConDepRatMatAbs::getMatrixType()
{
  return matrixType;
}



ostream &operator<<(ostream &os,const ConDepRatMatAbs &Q)
{
  Q.printOn(os);
  return os;
}



void ConDepRatMatAbs::printOn(ostream &os) const
{
  unsigned n=matrices.size();
  Sequence context;
  for(unsigned i=0 ; i<n ; ++i) {
    context.fromInt(i,order,alphabetMap);
    os<<context<<endl;
    RateMatrix *Q=matrices[i];
    if(!Q) continue;
    os<<*Q<<endl;
  }
  if(lowerOrderModel)
    lowerOrderModel->printOn(os);
}



int ConDepRatMatAbs::getNumMatrices() const {
  return matrices.size();
}



RateMatrix &ConDepRatMatAbs::getIthMatrix(int i) {
  return *matrices[i];
}



ConDepRatMatAbs *ConDepRatMatAbs::getLowerOrderModel(int order) {
    if(order<0 || order==this->order-1) return lowerOrderModel;
    return lowerOrderModel->getLowerOrderModel(order);
}



bool ConDepRatMatAbs::isDual() {
    return dual;
}



/*
  The following method changes the lower order models to make them
  consistent with the higher-order models.  It does this by averaging
  the higher-order RateMatrix's having the same lower-order context suffix.
 */
void ConDepRatMatAbs::averageLowerOrderModel(bool dual) {
    if(!lowerOrderModel) return;
    if(dual) {averageDual(); return;}
    Array1D<RateMatrix*> &lowerOrderMatrices=lowerOrderModel->matrices;
    int alphaRange=alphabetMap.getRangeSize();
    int numSuffixes=lowerOrderMatrices.size();
    for(int i=0 ; i<numSuffixes ; ++i) {
        Sequence suffix;
        suffix.fromInt(i,order-1,alphabetMap);
        Array1D<RateMatrix*> hiMs(alphaRange);
        for(Symbol y=0 ; y<alphaRange ; ++y) {
            Symbol x=alphabetMap.unmap(y);
            Sequence context;
            context.append(x);
            context.append(suffix);
            hiMs[y]=&lookup(context);
        }
        RateMatrix *&lowerOrderMatrix=lowerOrderMatrices[i];
        delete lowerOrderMatrix;
        lowerOrderMatrix=hiMs[0]->average(hiMs);
    }
    lowerOrderModel->averageLowerOrderModel(false);
}



void ConDepRatMatAbs::deleteOddOrders() {
    if(!lowerOrderModel) return;
    if(lowerOrderModel->order % 2) {
        ConDepRatMatAbs *newLowerOrderModel=lowerOrderModel->lowerOrderModel;
        lowerOrderModel->lowerOrderModel=NULL;
        delete lowerOrderModel;
        lowerOrderModel=newLowerOrderModel;
    }
    lowerOrderModel->deleteOddOrders();
}



void ConDepRatMatAbs::averageDual() {
    if(!lowerOrderModel) return;
    Array1D<RateMatrix*> &lowerOrderMatrices=lowerOrderModel->matrices;
    int halfOrder=order/2;
    int alphaRange=alphabetMap.getRangeSize();
    int numSuffixes=lowerOrderMatrices.size();
    Sequence parentSuffix, childSuffix, parentContext, childContext;
    for(int i=0 ; i<numSuffixes ; ++i) {
        parentSuffix.fromInt(i,halfOrder-1,alphabetMap);
        for(int j=0 ; j<numSuffixes ; ++j) {
            childSuffix.fromInt(j,halfOrder-1,alphabetMap);
            Array1D<RateMatrix*> hiMs(alphaRange*alphaRange);
            int nextM=0;
            for(Symbol w=0 ; w<alphaRange ; ++w) {
                Symbol x=alphabetMap.unmap(w);
                parentContext=x; parentContext.append(parentSuffix);
                for(Symbol y=0 ; y<alphaRange ; ++y) {
                    Symbol z=alphabetMap.unmap(y);
                    childContext=z; childContext.append(childSuffix);
                    hiMs[nextM++]=&lookup(parentContext,childContext);
                }
            }
            int lowerOrderMatrixIndex=
                lowerOrderModel->lookupIndex(parentSuffix,childSuffix);
            delete lowerOrderMatrices[lowerOrderMatrixIndex];
            lowerOrderMatrices[lowerOrderMatrixIndex]=hiMs[0]->average(hiMs);
        }
    }
    lowerOrderModel->averageDual();    
}



ConDepRatMatAbs *ConDepRatMatAbs::ancestralContextsOnly() {
    int acoOrder=order/2;
    ConDepRatMatAbs *newM=
        new ConDepRatMatAbs(acoOrder,seqType,matrixType,false);
    ancestralContextsOnly(*newM);
    return newM;
}


void ConDepRatMatAbs::ancestralContextsOnly(ConDepRatMatAbs &newM) {
    if(order==0) { // ZEROTH ORDER => JUST COPY
        newM.matrices[0]=matrices[0]->clone();
    }
    else if(order%2==0) { // EVEN ORDER => RECURSE TO LOWER ORDER
        lowerOrderModel->ancestralContextsOnly(newM);
    }
    else { // ODD ORDER => JUST COPY, THEN RECURSE
        int n=newM.matrices.size();
        for(int i=0 ; i<n ; ++i) {
            RateMatrix *Q=matrices[i];
            if(Q) newM.matrices[i]=Q->clone();
        }
        lowerOrderModel->ancestralContextsOnly(*newM.lowerOrderModel);
    }
}



ConDepRatMatAbs *ConDepRatMatAbs::clone() {
    ConDepRatMatAbs *M=new ConDepRatMatAbs(order,seqType,matrixType,dual);
    clone(*M);
    return M;
}



void ConDepRatMatAbs::clone(ConDepRatMatAbs &M) {
    int n=matrices.size();
    for(int i=0 ; i<n ; ++i) {
        RateMatrix *Q=matrices[i];
        M.matrices[i]=Q ? Q->clone() : NULL;
    }
    if(lowerOrderModel)
        lowerOrderModel->clone(*M.lowerOrderModel);
}



ConDepRatMatAbs *ConDepRatMatAbs::nextHigherOrder_dual() {
  // First, initialize the lower-order model within the new model, by
  // simply cloning this matrix into it
  int newOrder=order+2;
  ConDepRatMatAbs *M=new ConDepRatMatAbs(newOrder,seqType,matrixType,dual);
  ConDepRatMatAbs *Mlo=M->getLowerOrderModel(order);
  clone(*Mlo);

  // Now iterate through all contexts of the new model, and initialize each
  // by looking at a suffix (i.e., lower order)
  Sequence nmer;
  int numMatrices=M->matrices.size();
  for(int i=0 ; i<numMatrices ; ++i) {
    nmer.fromInt(i,newOrder,alphabetMap);
    int lowerOrderIndex=nmer.asInt(alphabetMap,2,order);
    RateMatrix *Q=matrices[lowerOrderIndex];
    M->matrices[i]=Q ? Q->clone() : NULL;
  }

  return M;
}



/********************************************************************
 This method returns a new model M of order N+1 in which all (N+1)-mers
 have their rate matrices initialized from the corresponding N-mer suffix
 within the current Nth-order model.  That is, this Nth-order model acts
 as a template for an (N+1)th order model, which will actually behave 
 identically to the lower order model (until perturbed, such as during
 training).
 */
ConDepRatMatAbs *ConDepRatMatAbs::nextHigherOrder() {
  if(dual) return nextHigherOrder_dual();

  // First, initialize the lower-order model within the new model, by
  // simply cloning this matrix into it
  int newOrder=order+1;
  ConDepRatMatAbs *M=new ConDepRatMatAbs(newOrder,seqType,matrixType,dual);
  ConDepRatMatAbs *Mlo=M->getLowerOrderModel(order);
  clone(*Mlo);

  // Now iterate through all contexts of the new model, and initialize each
  // by looking at a suffix (i.e., lower order)
  Sequence nmer;
  int numMatrices=M->matrices.size();
  for(int i=0 ; i<numMatrices ; ++i) {
    nmer.fromInt(i,newOrder,alphabetMap);
    int lowerOrderIndex=nmer.asInt(alphabetMap,1,order);
    M->matrices[i]=matrices[lowerOrderIndex]->clone();
  }

  return M;
}





