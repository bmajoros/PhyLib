/****************************************************************
 RateMatrixType.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_RateMatrixType_H
#define INCL_RateMatrixType_H
#include <iostream>
#include "BOOM/String.H"

enum RateMatrixType
  {
    MT_JC,   // Jukes-Cantor
    MT_FEL,  // Felsenstein
    MT_KIM,  // Kimura 2-parameter
    MT_HKY,  // HKY
    MT_REV,  // General reversible
    MT_NMER, // General reversible model for N-mers, assuming cond. indep.
    MT_AMINO,
    MT_CODON,
    MT_GAP,  // for the alphabet {N,-,.} ("N", "dash", "dot")
    MT_HB,   // Halpern-Bruno
    MT_GAINLOSS, // 2x2 matrix for gain/loss model
    MT_UNKNOWN
  };

ostream &operator<<(ostream &,RateMatrixType);
istream &operator>>(istream &,RateMatrixType &);
RateMatrixType rateMatrixTypeFromString(const BOOM::String &);

#endif
