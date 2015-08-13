/****************************************************************
 ContextType.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "ContextType.H"
using namespace std;
using namespace BOOM;

ostream &operator<<(ostream &os,ContextType t) 
{
  switch(t) 
    {
    case CT_HOG:   os<<"HOG"; break;
    case CT_ACO:   os<<"ACO"; break;
    case CT_RCO:   os<<"RCO"; break;
    case CT_LCO:   os<<"LCO"; break;
    case CT_TRCO:  os<<"TRCO"; break;
    case CT_MP:    os<<"MP"; break;
    case CT_NMER:  os<<"NMER"; break;
    case CT_AMINO: os<<"AMINO"; break;
    case CT_CODON: os<<"CODON"; break;
    default: throw "operator<<(os,ContextType)";
    }
  return os;
}



istream &operator>>(istream &is,ContextType &t) {
    String s;
    is>>s;
    s.trimWhitespace();
    t=contextTypeFromString(s);
    return is;
}



ContextType contextTypeFromString(const String &s) 
{
  if(s=="HOG") return CT_HOG;
  if(s=="ACO") return CT_ACO;
  if(s=="RCO") return CT_RCO;
  if(s=="LCO") return CT_LCO;
  if(s=="TRCO") return CT_TRCO;
  if(s=="MP") return CT_MP;
  if(s=="NMER") return CT_NMER;
  if(s=="AMINO") return CT_AMINO;
  if(s=="CODON") return CT_CODON;
  throw String("Unknown ContextType: ")+s;
}



