/****************************************************************
 predict-phase.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/Exceptions.H"
#include "BOOM/CommandLine.H"
#include "BOOM/MultiAlignment.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/Constants.H"
#include "Phylogeny.H"
#include "NthOrdRateMatrix.H"
#include "RateMatrixType.H"
#include "TRCO_Felsenstein.H"
#include "LCO_Felsenstein.H"
#include "RCO_Felsenstein.H"
#include "ACO_Felsenstein.H"
#include "HOG_Felsenstein.H"
#include "NmerFelsenstein.H"
#include "FitchFelsenstein.H"
#include "FitchParsimony.H"
#include "ContextType.H"
using namespace std;
using namespace BOOM;


class Application
{
  Alphabet *alphabet;
  BitSet gapSymbols;
  Phylogeny *phylogenies[3];
  NthOrdRateMatrix *rateMatrices[3];
  MultSeqAlignment *alignment;
  int alignmentPhase, alignmentLength;
  ContextType contextType;
  int numThreads;

  void loadModel(const String &filestem,int modelNum);
  void loadAlignment(const String &mafFile);
  double evaluate(Phylogeny *,NthOrdRateMatrix *,int alignmentPhase,
		  int modelPhase);
    double evaluate(int alignmentPhase);
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const String &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const RootException &e) 
      {
	cerr<<"BOOM exception caught in main: "<<e.getMessage()<<endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
  : contextType(CT_TRCO), alphabet(NULL), alignment(NULL), numThreads(1),
    gapSymbols(256)
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"c:");
    if(cmd.numArgs()!=4)
      throw String("predict-phase -c TRCO <model-filestem0> <model-filestem1> <model-filestem2> <*.maf>");
    String filestem[3];
    for(int i=0 ; i<3 ; ++i) filestem[i]=cmd.arg(i);
    String mafFilename=cmd.arg(3);
    if(cmd.option('c')) contextType=contextTypeFromString(cmd.optParm('c'));

    // Load models
    for(int i=0 ; i<3 ; ++i) loadModel(filestem[i],i);

    // Load the alignment
    alphabet=&DnaDashDotAlphabet::global();
    gapSymbols.addMember(alphabet->lookup('-'));
    gapSymbols.addMember(alphabet->lookup('.'));
    gapSymbols.addMember(alphabet->lookup('N'));
    loadAlignment(mafFilename);

    // Compute likelihoods and select phase
    int bestPhase=0;
    double bestLL=NEGATIVE_INFINITY;
    for(int phase=0 ; phase<3 ; ++phase) {
      double LL=evaluate(phase);
      if(LL>bestLL) {
	bestLL=LL;
	bestPhase=phase;
      }
    }
    cout<<bestPhase<<endl;

    return 0;
  }


void Application::loadModel(const String &filestem,int modelNum)
{
  String phyFile=filestem+".phy", matrixFile=filestem+".matrix";
  phylogenies[modelNum]=new Phylogeny(phyFile);
  rateMatrices[modelNum]=NthOrdRateMatrix::load(matrixFile);
}



void Application::loadAlignment(const String &mafFile)
{
  ifstream is(mafFile.c_str());
  if(!is.good()) throw String("Can't open file ")+mafFile;
  MultiAlignment *a=new MultiAlignment;
  a->loadMAF(is);
  is.close();
  if(a->getNumTracks()<2) {cout<<"no informants"<<endl;exit(0);}
  alignmentPhase=a->getPhase();
  alignmentLength=a->getLength();
  a->toupper();
  alignment=new MultSeqAlignment(*a,*alphabet,gapSymbols);
  delete a;
  for(int i=0 ; i<3 ; ++i) phylogenies[i]->attachAlignment(*alignment);
  alignment->deleteTargetGaps(phylogenies[0]->getRoot()->getID());
  if(contextType==CT_MP) {
    FitchParsimony fitch(*phylogenies[0],*alphabet,*alignment,
			 alignment->getGapSymbols());
    fitch.run();
  }
  for(int i=0 ; i<3 ; ++i) phylogenies[i]->attachAlignment(*alignment);
}



double Application::evaluate(int alignmentPhase)
{
  double LL=0;
  for(int modelPhase=0 ; modelPhase<3 ; ++modelPhase) {
    Phylogeny *phylogeny=phylogenies[modelPhase];
    NthOrdRateMatrix *Q=rateMatrices[modelPhase];
    LL+=evaluate(phylogeny,Q,alignmentPhase,modelPhase);
  }
  return LL;
}



double Application::evaluate(Phylogeny *phylogeny,NthOrdRateMatrix *Q,
			     int alignmentPhase,int modelPhase)
{
  double LL=NEGATIVE_INFINITY;
  int order=Q->getOrder();
  bool dual=Q->isDual();
  int phase=(modelPhase+3-alignmentPhase)%3;
  switch(contextType) 
    {
    case CT_MP: {
      if(dual) {
	Q=Q->ancestralContextsOnly();
	order/=2;
      }
      FitchFelsenstein fel(order,*phylogeny,*Q,*alignment,DNA,numThreads);
      LL=fel.logLikelihoodInPhase(0,alignmentLength,phase);
    }
      break;
    case CT_LCO: {
      if(dual) {
        Q=Q->ancestralContextsOnly();
        order/=2;
      }
      LCO_Felsenstein fel(order,*phylogeny,*Q,*alignment,DNA,numThreads);
      LL=fel.logLikelihoodInPhase(0,alignmentLength,phase);
    }
      break;
    case CT_TRCO: {
      if(dual) {
	Q=Q->ancestralContextsOnly();
	order/=2;
      }
      TRCO_Felsenstein fel(order,*phylogeny,*Q,*alignment,DNA,numThreads);
      LL=fel.logLikelihoodInPhase(0,alignmentLength,phase);
    }
      break;
    case CT_RCO: {
      if(dual) {
	Q=Q->ancestralContextsOnly();
	order/=2;
      }
      RCO_Felsenstein fel(order,*phylogeny,*Q,*alignment,DNA,numThreads);
      LL=fel.logLikelihoodInPhase(0,alignmentLength,phase);
    }
      break;
    case CT_ACO: {
      if(dual) {
	Q=Q->ancestralContextsOnly();
	order/=2;
      }
      ACO_Felsenstein fel(order,*phylogeny,*Q,*alignment,DNA,numThreads);
      LL=fel.logLikelihoodInPhase(0,alignmentLength,phase);
    }
      break;
    case CT_HOG: {
      throw "not implemented";
      //if(!dual) throw "HOG model requires dual contexts";
      //HOG_Felsenstein fel(order/2,*phylogeny,*Q,*alignment,DNA,nmerTable,
      //		  numThreads);
      //LL=fel.logLikelihoodInPhase(0,alignmentLength,phase);
    }
      break;
    case CT_NMER: {
      throw "not implemented";
      //NmerRateMatrix *nmerQ=static_cast<NmerRateMatrix*>(Q);
      //NmerFelsenstein fel(*phylogeny,nmerQ,*alignment,nmerTable,numThreads);
      //LL=fel.logLikelihoodInPhase(0,alignmentLength,phase);
      break;
    }
    default: throw "unrecognized context type";
    }
  return LL;
}
