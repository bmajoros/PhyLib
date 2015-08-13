/****************************************************************
 RateMatrixType.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/

#include "RateMatrixType.H"


ostream &operator<<(ostream &os,RateMatrixType t)
{
  switch(t)
    {
    case MT_JC:   os<<"JC";  break;
    case MT_FEL:  os<<"FEL"; break;
    case MT_KIM:  os<<"KIM"; break;
    case MT_HKY:  os<<"HKY"; break;
    case MT_REV:  os<<"REV"; break;
    case MT_NMER: os<<"NMER"; break;
    case MT_AMINO:os<<"AMINO"; break;
    case MT_CODON:os<<"CODON"; break;
    case MT_GAP:  os<<"GAP";   break;
    case MT_HB:   os<<"HB";  break;
    case MT_GAINLOSS:os<<"GAINLOSS"; break;
    }
  return os;
}



RateMatrixType rateMatrixTypeFromString(const BOOM::String &s)
{
  if(s=="JC")  return MT_JC;
  if(s=="FEL") return MT_FEL;
  if(s=="KIM") return MT_KIM;
  if(s=="HKY") return MT_HKY;
  if(s=="REV") return MT_REV;
  if(s=="NMER") return MT_NMER;
  if(s=="AMINO") return MT_AMINO;
  if(s=="CODON") return MT_CODON;
  if(s=="GAP") return MT_GAP;
  if(s=="HB")  return MT_HB;
  if(s=="GAINLOSS") return MT_GAINLOSS;
  throw BOOM::String("Unknown rate matrix type: ")+s;
}



istream &operator>>(istream &is,RateMatrixType &t)
{
  BOOM::String s;
  is>>s;
  t=rateMatrixTypeFromString(s);
  return is;
}


