/****************************************************************
 AminoRateMatrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "AminoRateMatrix.H"
using namespace std;
using namespace BOOM;

AminoRateMatrix::AminoRateMatrix()
    : RateMatrix(PROTEIN,MT_AMINO)
{
  // ctor
}



void AminoRateMatrix::save(ostream &os) {
}



RateMatrixType AminoRateMatrix::getType() {
}



int AminoRateMatrix::numParameters() const {
}



double AminoRateMatrix::getIthParm(int i) const {
}



void AminoRateMatrix::setIthParm(int i,double parm) {
}



void AminoRateMatrix::partialDerivative(int parameter,RateMatrix &) const {
}



RateMatrix *AminoRateMatrix::average(const Array1D<RateMatrix*> &) const {
}



void AminoRateMatrix::addNoise(GSL::ContinuousDistribution &) {
}



RateMatrix *AminoRateMatrix::clone() const {
}







