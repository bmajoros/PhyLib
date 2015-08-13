/****************************************************************
 GapRewritingModel.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "GapRewritingModel.H"
#include "RCO_Felsenstein.H"
#include "LCO_Felsenstein.H"

using namespace std;
using namespace BOOM;


GapRewritingModel::GapRewritingModel(int order,
				     Phylogeny &phylogeny,
				     NthOrdRateMatrix &gapModel,
				     NthOrdRateMatrix &nucModel)
  : gapModel(gapModel),
    nucModel(nucModel),
    phylogeny(phylogeny),
    order(order)
{
}



double GapRewritingModel::logLikelihood(MultSeqAlignment &alignment,
					int column)
{
  RCO_Felsenstein f1(order,phylogeny,nucModel,alignment,DNA);
  LCO_Felsenstein f2(order,phylogeny,gapModel,alignment,DNA);
  return f1.logLikelihood(column)+f2.logLikelihood(column);
}



double GapRewritingModel::logLikelihood(MultSeqAlignment &alignment)
{
  RCO_Felsenstein f1(order,phylogeny,nucModel,alignment,DNA);
  LCO_Felsenstein f2(order,phylogeny,gapModel,alignment,DNA);
  return f1.logLikelihood()+f2.logLikelihood();
}



double GapRewritingModel::logLikelihood(MultSeqAlignment &alignment,
					int begin,int end)
{
  RCO_Felsenstein f1(order,phylogeny,nucModel,alignment,DNA);
  LCO_Felsenstein f2(order,phylogeny,gapModel,alignment,DNA);
  return f1.logLikelihood(begin,end)+f2.logLikelihood(begin,end);
}


