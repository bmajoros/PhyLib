/****************************************************************
 RateMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_RateMatrix_H
#define INCL_RateMatrix_H
#include <iostream>
#include "BOOM/GSL/Matrix.H"
#include "BOOM/MolecularSequenceType.H"
#include "BOOM/AlphabetMap.H"
#include "BOOM/GSL/ContinuousDistribution.H"
#include "SubstitutionMatrix.H"
#include "RateMatrixType.H"
using namespace std;
using namespace BOOM;


void rateMatrixToSubstMatrix(const GSL::Matrix &rateMatrix,
			     double branchLength,
			     GSL::Matrix &substMatrix);


/****************************************************************
 RateMatrix : a matrix of instantaneous substitution rates, usu.
 denoted in the literature as Q.  Given a branch length t, this
 rate matrix can be instantiated into a substitution matrix, P(t).

 Note that for DNA matrices, only {A,C,G,T} are represented, so
 that the 'N' character is not allocated a row/column in the matrix.
 ****************************************************************/
class RateMatrix
{
public:
  RateMatrix(MolecularSequenceType,RateMatrixType=MT_UNKNOWN);
  RateMatrix(MolecularSequenceType,Alphabet &,RateMatrixType=MT_UNKNOWN);
  virtual ~RateMatrix() {}
  static RateMatrix *load(const BOOM::String &filename);
  static RateMatrix *load(istream &);
  static RateMatrix *random(MolecularSequenceType,RateMatrixType,
			    BOOM::Array1D<double> equilibriumFreqs=omitted);
  void save(const BOOM::String &filename);
  virtual void save(ostream &)=0;
  virtual RateMatrixType getType() {return matrixType;}
  double operator()(Symbol row,Symbol column) const;
  double &operator()(Symbol row,Symbol column);
  SubstitutionMatrix *instantiate(double branchLength);
  void instantiate(double branchLength,SubstitutionMatrix &);
  int getNumRows() const;
  const Alphabet &getAlphabet() const;
  AlphabetMap &getAlphabetMap() {return alphabetMap;}
  virtual int numParameters() const;
  virtual double getIthParm(int i) const=0;
  virtual void setIthParm(int i,double parm)=0;
  virtual void partialDerivative(int parameter,RateMatrix &) const;
  void initMatrix_ACGT(double *array);
  void installDiagonals(); // makes sure rows sum to zero
  RateMatrixType getMatrixType();
  virtual RateMatrix *average(const Array1D<RateMatrix*> &) const=0;
  virtual void addNoise(GSL::ContinuousDistribution &)=0;
  virtual RateMatrix *clone() const=0;
  virtual void init()=0;
  void rescale(double factor);
  const GSL::Matrix &peek() const;
  RateMatrix *complement() const;
protected:
  RateMatrixType matrixType;
  static BOOM::Array1D<double> omitted;
  Alphabet &alphabet;
  GSL::Matrix M;
  MolecularSequenceType seqType;
  AlphabetMap alphabetMap;
  
  friend ostream &operator<<(ostream &,const RateMatrix &);
  void printOn(ostream &) const;
};

ostream &operator<<(ostream &,const SubstitutionMatrix &);


#endif
