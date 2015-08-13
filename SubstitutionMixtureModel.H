/****************************************************************
 SubstitutionMixtureModel.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_SubstitutionMixtureModel_H
#define INCL_SubstitutionMixtureModel_H
#include <iostream>
#include "RateMatrixDecomposed.H"
using namespace std;


/****************************************************************
 SubstitutionMixtureModel = a mixture between two rate matrices
 along a phylogeny branch. Let t be the branch length and s a
 part of that branch (s<t). If the ancestral rate matrix was Q1,
 but that matrix changed at some point to Q2, the expectation of
 the resulting matrix is the integral of P1(s)*P2(t-s)ds between
 0 and t.  That's what this class computes.

 Notes: (1) reversing Q1 and Q2 does change the result. (2) this
 implementation is quite fast, and several heuristics were found
 to perform less well---such as simply averaging P1(t) and P2(t),
 or simplying taking P1(t/2)*P2(t/2).
 ****************************************************************/
struct SubstitutionMixtureModel : public SubstitutionMatrix {
  SubstitutionMixtureModel(const RateMatrix &Q1,const RateMatrix &Q2,
			   double branchLength);
  SubstitutionMixtureModel(const RateMatrixDecomposed &Q1,
			   const RateMatrixDecomposed &Q2,
			   double branchLength);
protected:
  void init(const RateMatrixDecomposed &Q1,
	    const RateMatrixDecomposed &Q2,
	    double branchLength);
};

#endif

