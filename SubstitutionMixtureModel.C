/****************************************************************
 SubstitutionMixtureModel.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "SubstitutionMixtureModel.H"
using namespace std;
using namespace BOOM;

const double TOLERANCE=0.0001;


SubstitutionMixtureModel::SubstitutionMixtureModel(const RateMatrix &Q1,
						   const RateMatrix &Q2,
						   double t)
  : SubstitutionMatrix(Q1.getAlphabet(),Q1.getAlphabetMap())
{
  init(RateMatrixDecomposed(Q1),RateMatrixDecomposed(Q2),t);
}



SubstitutionMixtureModel::SubstitutionMixtureModel(const 
                RateMatrixDecomposed &Q1,const RateMatrixDecomposed &Q2,
						   double t)
  : SubstitutionMatrix(Q1.getAlphabet(),Q1.getAlphabetMap())
{
  init(Q1,Q2,t);
}



void SubstitutionMixtureModel::init(const RateMatrixDecomposed &Q1,
				    const RateMatrixDecomposed &Q2,
				    double t)
{
  const GSL::Matrix &M1=Q1.getLeftMatrix(), &M1inv=Q1.getRightMatrix();
  const GSL::Matrix &M2=Q2.getLeftMatrix(), &M2inv=Q2.getRightMatrix();
  const GSL::Vector &lambda1=Q1.getEigenvalues();
  const GSL::Vector &lambda2=Q2.getEigenvalues();
  int n=M1.getNumRows();

  // Form the "m" matrix, m = M1inv * M2
  GSL::Matrix m;
  M1inv.times(M2,m);

  // Form the middle matrix, called "alpha"
  GSL::Matrix alpha(n,n);
  for(int i=0 ; i<n ; ++i) {
    double expLambda2i_x_t=exp(lambda2[i]*t);
    for(int j=0 ; j<n ; ++j) {
      double diff=lambda1[j]-lambda2[i];
      if(fabs(diff)<TOLERANCE)
	alpha(j,i)=m(j,i)*expLambda2i_x_t;
      else
	alpha(j,i)=m(j,i)*(exp(lambda1[j]*t)-expLambda2i_x_t)/(diff*t);
     }
  }

  // Multiple by M1 on the left
  GSL::Matrix temp;
  M1.times(alpha,temp);

  // Multiply by M2inv on the right
  temp.times(M2inv,M);
}



