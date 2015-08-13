/****************************************************************
 ContextType.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ContextType_H
#define INCL_ContextType_H
#include <iostream>
#include "BOOM/String.H"
using namespace std;
using namespace BOOM;


enum ContextType {
    CT_HOG,    // full context modeling in both ancestor and descendent
    CT_ACO,    // Ancestral Contexts Only
    CT_RCO,    // Root Contexts Only
    CT_TRCO,   // Transitive Root Contexts Only
    CT_LCO,    // Leaf Contexts Only
    CT_MP,     // Maximum Parsimony
    CT_NMER,   // N-mer substitution model (1 change per generation only)
    CT_AMINO,  // Amino Acid substitution model
    CT_CODON   // Codon substitution model (using AMINO model)
};

ostream &operator<<(ostream &,ContextType);
istream &operator>>(istream &,ContextType &);
ContextType contextTypeFromString(const String &);

#endif

