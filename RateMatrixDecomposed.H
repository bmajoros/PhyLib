/****************************************************************
 RateMatrixDecomposed.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_RateMatrixDecomposed_H
#define INCL_RateMatrixDecomposed_H
#include <iostream>
#include "RateMatrix.H"
#include "BOOM/GSL/Vector.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 RateMatrixDecomposed : a rate matrix that has been decomposed 
 via Q=ABC, for A the matrix of column eigenvectors, C the inverse
 of A, and B the diagonal matrix of eigenvalues of Q.  Computing
 the substitution matrix P(t) for many t's is faster using this
 class than using RateMatrix, since RateMatrix will repeat the
 (identical) spectral decomposition for each t, which is unnecessary.
 ****************************************************************/
struct RateMatrixDecomposed {
  RateMatrixDecomposed(const RateMatrix &Q);
  RateMatrixDecomposed(const GSL::Matrix &left,const GSL::Vector &eigenvalues,
		       const GSL::Matrix &right,const Alphabet &,
		       const AlphabetMap &);
  const GSL::Matrix &getLeftMatrix() const;
  const GSL::Matrix &getRightMatrix() const;
  const GSL::Vector &getEigenvalues() const;
  void getEigenvalueMatrix(GSL::Matrix &) const; // diag(lambda)
  void getExpEigenMatrix(GSL::Matrix &,double t) const;// diag(exp(lambda*t))
  const Alphabet &getAlphabet() const;
  const AlphabetMap &getAlphabetMap() const;
  SubstitutionMatrix *instantiate(double t) const;
protected:
  GSL::Matrix left, right; // Q = left * diag(eigenvalues) * right
  GSL::Vector eigenvalues;
  const Alphabet &alphabet;
  const AlphabetMap &alphabetMap;
  void init(const RateMatrix &);
};

#endif

