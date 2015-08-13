/****************************************************************
 RateMatrixDecomposed.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "RateMatrixDecomposed.H"
using namespace std;
using namespace BOOM;


RateMatrixDecomposed::RateMatrixDecomposed(const RateMatrix &Q)
  : alphabet(Q.getAlphabet()), alphabetMap(Q.getAlphabetMap())
{
  init(Q);
}



RateMatrixDecomposed::RateMatrixDecomposed(const GSL::Matrix &left,
					   const GSL::Vector &eigenvalues,
					   const GSL::Matrix &right,
					   const Alphabet &alphabet,
					   const AlphabetMap &alphabetMap)
  : left(left), eigenvalues(eigenvalues), right(right),
    alphabet(alphabet), alphabetMap(alphabetMap)
{
  // ctor
}



const GSL::Matrix &RateMatrixDecomposed::getLeftMatrix() const
{
  return left;
}



const GSL::Matrix &RateMatrixDecomposed::getRightMatrix() const
{
  return right;
}



const GSL::Vector &RateMatrixDecomposed::getEigenvalues() const
{
  return eigenvalues;
}



void RateMatrixDecomposed::init(const RateMatrix &Q)
{
  GSL::Vector imagParts; // should be all zeros for practical rate matrices
  Q.peek().getEigenVectors(left,eigenvalues,imagParts);
  left.invert(right);
}



void RateMatrixDecomposed::getEigenvalueMatrix(GSL::Matrix &M) const
{
  int n=left.getNumRows();
  M.resize(n,n);
  M.setAllTo(0.0);
  for(int i=0 ; i<n ; ++i)
    M(i,i)=exp(eigenvalues[i]);
}



void RateMatrixDecomposed::getExpEigenMatrix(GSL::Matrix &M,double t) const
{
  int n=left.getNumRows();
  M.resize(n,n);
  M.setAllTo(0.0);
  for(int i=0 ; i<n ; ++i)
    M(i,i)=exp(eigenvalues[i]*t);
}



const Alphabet &RateMatrixDecomposed::getAlphabet() const
{
  return alphabet;
}



const AlphabetMap &RateMatrixDecomposed::getAlphabetMap() const
{
  return alphabetMap;
}



SubstitutionMatrix *RateMatrixDecomposed::instantiate(double t) const
{
  GSL::Matrix lambda, leftLambda;
  getExpEigenMatrix(lambda,t);
  left.times(lambda,leftLambda);
  SubstitutionMatrix *Pt=new SubstitutionMatrix(alphabet,alphabetMap);
  leftLambda.times(right,Pt->peek());
  return Pt;
}

