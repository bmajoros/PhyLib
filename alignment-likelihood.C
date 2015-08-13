
/****************************************************************
 alignment-likelihood.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/MultiAlignment.H"
#include "BOOM/Vector.H"
#include "BOOM/Random.H"
#include "BOOM/Constants.H"
#include "BOOM/Environment.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/Time.H"
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
#include "AlignmentNmerTable.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                       class Application
 ****************************************************************/
class Application
{
  Alphabet *alphabet;
  AlignmentNmerTable *nmerTable;
  bool wobble, bestFrameOnly, onePhaseOnly;
  int desiredPhase;
  Phylogeny *phylogeny;
  MultSeqAlignment *alignment;
  BitSet gapSymbols;
  ContextType contextType;
  void loadAlignments(const String &mafFile,bool shouldSlice);
public:
  Application();
  int main(int argc,char *argv[]);
};



/****************************************************************
                              main()
 ****************************************************************/
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
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
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




/****************************************************************
                      Application methods
 ****************************************************************/



Application::Application()
  : gapSymbols(256)
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"c:p:o:wfz:");
  if(cmd.numArgs()!=3)
    throw String(
"\nalignment-likelihood [-cpow] <*.phy> <*.maf> <*.matrix>\n\
\n\
   where -c <context-type> = HOG | RCO | TRCO | MP | AMINO | CODON\n\
         -p <N> = evaluate only a prefix of length <N>\n\
         -o <N> = use order N (must be less than model's actual order)\n\
         -w = evaluate only the 3rd codon position (\"wobble\")\n\
         -f = evaluate in all three frames and report only the highest\n\
         -z <N> = evaluate only phase N\n\
\n\n");
  String phyFile=cmd.arg(0);
  String mafFile=cmd.arg(1);
  String matrixFile=cmd.arg(2);
  int numThreads=Environment::lookup("NUM_CPUS");
  if(numThreads<1) numThreads=1;
  wobble=cmd.option('w');
  bestFrameOnly=cmd.option('f');
  onePhaseOnly=cmd.option('z');
  if(onePhaseOnly) desiredPhase=cmd.optParm('z').asInt();

  // Load the phylogeny
  phylogeny=new Phylogeny(phyFile);
  //phylogeny->attachAlignment(*alignment);

  // Load the alignments
  alphabet=&DnaDashDotAlphabet::global();
  gapSymbols.addMember(alphabet->lookup('-'));
  gapSymbols.addMember(alphabet->lookup('.'));
  gapSymbols.addMember(alphabet->lookup('N'));
  loadAlignments(mafFile,onePhaseOnly || wobble);
  /*
  Vector<MultiAlignment*> alignments;
  MultiAlignment::loadMAF(mafFile,alignments);
  MultiAlignment *a=MultiAlignment::combine(alignments,true);
  int phase=a->getPhase();
  a->toupper();
  //gapSymbols.addMember(DnaAlphabet::global().lookup('N'));
  MultSeqAlignment MSA(*a,alphabet,gapSymbols);
  MultSeqAlignment *alignment=&MSA;
  if((wobble||onePhaseOnly) && phase>0) {
    alignment=MSA.getSlice(3-phase,MSA.getLength()-3+phase);
  }
  */

  // Load the rate matrix
  NthOrdRateMatrix *Q=NthOrdRateMatrix::load(matrixFile);
  int order=Q->getOrder();
  bool dual=Q->isDual();
  if(cmd.option('o')) {
    int O=cmd.optParm('o').asInt();
    if(O>=order) throw "option -o can only reduce order, not increase it";
    Q=Q->getLowerOrderModel(O);
    order=Q->getOrder();
  }

  // Delete alignment columns having a gap in the root sequence
  int rootID=phylogeny->getRoot()->getID();
  alignment->deleteTargetGaps(rootID);
  Q->getAlphabetMap();
  nmerTable=new AlignmentNmerTable(*alignment,Q->getAlphabetMap(),order);
  int length=alignment->getLength();
  if(cmd.option('p')) {
    int p=cmd.optParm('p').asInt();
    if(p<length) length=p;
  }

  // Identify the desired context type
  contextType=CT_TRCO;
  if(cmd.option('c')) contextType=contextTypeFromString(cmd.optParm('c'));

  // Compute likelihood
  Time timer;
  timer.startCounting();
  double LL=0;
  switch(contextType) 
    {
    case CT_MP: {
      if(dual) {
	Q=Q->ancestralContextsOnly();
	order/=2;
      }
      FitchParsimony fitch(*phylogeny,*alphabet,*alignment,gapSymbols);
      fitch.run();
      //phylogeny->attachAlignment(A);
      FitchFelsenstein fel(order,*phylogeny,*Q,*alignment,DNA,numThreads);
      if(wobble) LL=fel.logLikelihood3(0,length);
      else if(bestFrameOnly) LL=fel.logLikelihood_bestFrame(0,length);
      else if(onePhaseOnly) LL=fel.logLikelihoodInPhase(0,length,desiredPhase);
      else       LL=fel.logLikelihood(0,length);
    }
      break;

    case CT_LCO: {
      if(dual) {
        Q=Q->ancestralContextsOnly();
        order/=2;
      }
      LCO_Felsenstein fel(order,*phylogeny,*Q,*alignment,DNA,numThreads);
      if(wobble) LL=fel.logLikelihood3(0,length);
      else if(bestFrameOnly) LL=fel.logLikelihood_bestFrame(0,length);
      else if(onePhaseOnly) LL=fel.logLikelihoodInPhase(0,length,desiredPhase);
      else       LL=fel.logLikelihood(0,length);
    }
      break;
    case CT_TRCO: {
      if(dual) {
	Q=Q->ancestralContextsOnly();
	order/=2;
      }
      TRCO_Felsenstein fel(order,*phylogeny,*Q,*alignment,DNA,numThreads);
      if(wobble) LL=fel.logLikelihood3(0,length);
      else if(bestFrameOnly) LL=fel.logLikelihood_bestFrame(0,length);
      else if(onePhaseOnly) LL=fel.logLikelihoodInPhase(0,length,desiredPhase);
      else       LL=fel.logLikelihood(0,length);
    }
      break;
    case CT_RCO: {
      if(dual) {
	Q=Q->ancestralContextsOnly();
	order/=2;
      }
      RCO_Felsenstein fel(order,*phylogeny,*Q,*alignment,DNA,numThreads);
      if(wobble) LL=fel.logLikelihood3(0,length);
      else if(bestFrameOnly) LL=fel.logLikelihood_bestFrame(0,length);
      else if(onePhaseOnly) LL=fel.logLikelihoodInPhase(0,length,desiredPhase);
      else       LL=fel.logLikelihood(0,length);
    }
      break;
    case CT_ACO: {
      if(dual) {
	Q=Q->ancestralContextsOnly();
	order/=2;
      }
      //phylogeny->attachAlignment(A);
      ACO_Felsenstein fel(order,*phylogeny,*Q,*alignment,DNA,numThreads);
      if(wobble) LL=fel.logLikelihood3(0,length);
      else if(bestFrameOnly) LL=fel.logLikelihood_bestFrame(0,length);
      else if(onePhaseOnly) LL=fel.logLikelihoodInPhase(0,length,desiredPhase);
      else       LL=fel.logLikelihood(0,length);
    }
      break;
    case CT_HOG: {
      if(!dual) throw "HOG model requires dual contexts";
      //phylogeny->attachAlignment(A);
      HOG_Felsenstein fel(order/2,*phylogeny,*Q,*alignment,DNA,nmerTable,
			  numThreads);
      if(wobble) LL=fel.logLikelihood3(0,length);
      else if(bestFrameOnly) LL=fel.logLikelihood_bestFrame(0,length);
      else if(onePhaseOnly) LL=fel.logLikelihoodInPhase(0,length,desiredPhase);
      else       LL=fel.logLikelihood(0,length);
    }
      break;
    case CT_NMER: {
      NmerRateMatrix *nmerQ=static_cast<NmerRateMatrix*>(Q);
      NmerFelsenstein fel(*phylogeny,nmerQ,*alignment,nmerTable,numThreads);
      if(wobble) LL=fel.logLikelihood3(0,length);
      else if(bestFrameOnly) LL=fel.logLikelihood_bestFrame(0,length);
      else if(onePhaseOnly) LL=fel.logLikelihoodInPhase(0,length,desiredPhase);
      else       LL=fel.logLikelihood(0,length);
      break;
    }
    default: throw "unrecognized context type";
    }
  timer.stopCounting();
  
  // Generate output
  cout<<LL<<endl;
  cerr<<timer.elapsedTime()<<" "<<length<<" columns"<<endl;
  
  return 0;
}



void Application::loadAlignments(const String &mafFile,bool shouldSlice)
{
  alignment=NULL;
  ifstream is(mafFile.c_str());
  if(!is.good()) throw String("Can't open file ")+mafFile;
  while(!is.eof()) {
    MultiAlignment *a=new MultiAlignment;
    a->loadMAF(is);
    if(a->getNumTracks()>0) {
      a->toupper();
      MultSeqAlignment *m=new MultSeqAlignment(*a,*alphabet,gapSymbols);
      phylogeny->attachAlignment(*m);
      m->deleteTargetGaps(phylogeny->getRoot()->getID());
      MultSeqAlignment *slice=m;
      if(shouldSlice) {
	int phase=a->getPhase();
	int len=a->getLength();
	int begin=(3-phase)%3;
	int end=len;
	int endPhase=(phase+len)%3;
	end-=endPhase;
	int sliceLen=end-begin;
	if(sliceLen%3>0) throw "Application::loadAlignments";//###DEBUGGING
	slice=m->getSlice(begin,sliceLen);
	delete m;
      }
      if(!alignment) alignment=slice;
      else {
	alignment->append(*slice);
	delete slice;
      }
    }
    delete a;
  }
  is.close();
  if(contextType==CT_MP) {
    cout<<"Inferring ancestral sequences via parsimony..."<<endl;
    FitchParsimony fitch(*phylogeny,*alphabet,*alignment,
			 alignment->getGapSymbols());
    fitch.run();
  }
  phylogeny->attachAlignment(*alignment);
  alignment->setPhase(0);
}
