/***********************************************************************
 SubstitutionMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3
 (see www.opensource.org).
 ***********************************************************************/
#ifndef INCL_SubstitutionMatrix_H
#define INCL_SubstitutionMatrix_H
#include <iostream>
#include "BOOM/GSL/Matrix.H"
#include "BOOM/MolecularSequenceType.H"
#include "BOOM/Array1D.H"
#include "BOOM/AlphabetMap.H"
#include "BOOM/RouletteWheel.H"
using namespace BOOM;

/****************************************************************
 SubstitutionMatrix : a matrix giving the probability of given
 single-letter substitutions.

 Note that for DNA matrices, only {A,C,G,T} are represented, so
 that the 'N' character is not allocated a row/column in the matrix.
 ****************************************************************/
class SubstitutionMatrix
{
public:
  SubstitutionMatrix(const String &filename,MolecularSequenceType);
  SubstitutionMatrix(const Alphabet &,const AlphabetMap &);
  SubstitutionMatrix(const SubstitutionMatrix &);
  GSL::Matrix &peek();
  double operator()(Symbol row,Symbol column) const;
  double &operator()(Symbol row,Symbol column);
  void convertToLogs();
  bool isInLogSpace() const;
  void initGenerators();
  Symbol mutate(Symbol);
  void getEqFreqs(Array1D<double> &); // func resizes array for you
  const Alphabet &getAlphabet() {return alphabet;}
  const AlphabetMap &getAlphabetMap() {return alphabetMap;}
  void setAlphabetMap(const AlphabetMap &);
protected:
  const Alphabet &alphabet;
  AlphabetMap alphabetMap;
  GSL::Matrix M;
  MolecularSequenceType seqType;
  Array1D<RouletteWheel> generators; // for evolutionary simulations
  bool logSpace;
  
  friend ostream &operator<<(ostream &,const SubstitutionMatrix &);
  void printOn(ostream &) const;
};

ostream &operator<<(ostream &,const SubstitutionMatrix &);



#endif
