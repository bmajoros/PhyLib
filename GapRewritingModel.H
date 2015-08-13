/****************************************************************
 GapRewritingModel.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_GapRewritingModel_H
#define INCL_GapRewritingModel_H
#include <iostream>
#include "NthOrdRateMatrix.H"
#include "Phylogeny.H"
#include "BOOM/MultSeqAlignment.H"
using namespace std;
using namespace BOOM;


class GapRewritingModel
{
public:
  GapRewritingModel(int order,Phylogeny &,
		    NthOrdRateMatrix &gapModel,
		    NthOrdRateMatrix &nucModel);
  virtual double logLikelihood(MultSeqAlignment &,int column);
  virtual double logLikelihood(MultSeqAlignment &);
  virtual double logLikelihood(MultSeqAlignment &,int begin,int end);
protected:
  NthOrdRateMatrix &gapModel;
  NthOrdRateMatrix &nucModel;
  Phylogeny &phylogeny;
  int order;
};

#endif

