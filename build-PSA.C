/****************************************************************
 build-PSA.C
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
#include "BOOM/DnaAlphabet.H"
#include "BOOM/File.H"
#include "Phylogeny.H"
#include "NthOrdRateMatrix.H"
#include "RateMatrixType.H"
#include "TRCO_Felsenstein.H"
#include "RCO_Felsenstein.H"
#include "ACO_Felsenstein.H"
#include "HOG_Felsenstein.H"
#include "NmerFelsenstein.H"
#include "FitchFelsenstein.H"
#include "FitchParsimony.H"
#include "ContextType.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                       class Application
 ****************************************************************/
class Application
{
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
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"c:p:");
    if(cmd.numArgs()!=4)
        throw String(
"\nbuild-PSA [-cp] <*.phy> <*.maf> <*.matrix> <outfile>\n\
\n\
   where -c <context-type> = HOG | RCO | TRCO | MP | AMINO | CODON\n\
         -p <N> = evaluate only a prefix of length <N>\n\
\n\n");
    String phyFile=cmd.arg(0);
    String mafFile=cmd.arg(1);
    String matrixFile=cmd.arg(2);
    String outfile=cmd.arg(3);

    int numThreads=Environment::lookup("NUM_CPUS");
    if(numThreads<1) numThreads=1;
    
    // Load the alignments
    Vector<MultiAlignment*> alignments;
    MultiAlignment::loadMAF(mafFile,alignments);
    MultiAlignment *alignment=MultiAlignment::combine(alignments,true);
    alignment->toupper();
    Symbol gapSymbol=DnaAlphabet::global.lookup('N');
    MultSeqAlignment A(*alignment,DnaAlphabet::global,gapSymbol);
    int length=A.getLength();
    if(cmd.option('p')) {
      int p=cmd.optParm('p').asInt();
      if(p<length) length=p;
    }

    // Load the phylogeny
    Phylogeny *phylogeny=new Phylogeny(phyFile);
    phylogeny->attachAlignment(A);

    // Delete alignment columns having a gap in the root sequence
    int rootID=phylogeny->getRoot()->getID();
    A.deleteTargetGaps(rootID);

    // Load the rate matrix
    NthOrdRateMatrix *Q=NthOrdRateMatrix::load(matrixFile);
    int order=Q->getOrder();
    bool dual=Q->isDual();

    // Identify the desired context type
    ContextType contextType=CT_MP;
    if(cmd.option('c')) contextType=contextTypeFromString(cmd.optParm('c'));
    
    // Allocate appropriate version of Felsenstein's algorithm
    FelsensteinInterface *fel;
    switch(contextType) 
      {
      case CT_MP: {
	if(dual) {
	  Q=Q->ancestralContextsOnly();
	  order/=2;
	}
	FitchParsimony fitch(*phylogeny,DnaAlphabet::global,A,gapSymbol);
	fitch.run();
	phylogeny->attachAlignment(A);
	fel=new FitchFelsenstein(order,*phylogeny,*Q,A,DNA,numThreads);
      }
	break;
      case CT_TRCO: {
	if(dual) {
	  Q=Q->ancestralContextsOnly();
	  order/=2;
	}
	phylogeny->attachAlignment(A);
	fel=new TRCO_Felsenstein(order,*phylogeny,*Q,A,DNA,numThreads);
      }
	break;
      case CT_RCO: {
	if(dual) {
	  Q=Q->ancestralContextsOnly();
	  order/=2;
	}
	phylogeny->attachAlignment(A);
	fel=new RCO_Felsenstein(order,*phylogeny,*Q,A,DNA,numThreads);
      }
	break;
      case CT_ACO: {
	if(dual) {
	  Q=Q->ancestralContextsOnly();
	  order/=2;
	}
	phylogeny->attachAlignment(A);
	fel=new ACO_Felsenstein(order,*phylogeny,*Q,A,DNA,numThreads);
      }
	break;
      case CT_HOG: {
	if(!dual) throw "HOG model requires dual contexts";
	phylogeny->attachAlignment(A);
	fel=new HOG_Felsenstein(order/2,*phylogeny,*Q,A,DNA,numThreads);
      }
	break;
      case CT_NMER: {
	NmerRateMatrix *nmerQ=static_cast<NmerRateMatrix*>(Q);
	fel=new NmerFelsenstein(*phylogeny,nmerQ,A,numThreads);
	break;
      }
      default: throw "unrecognized context type";
      }

    // Write prefix-sum-array
    File file(outfile,"w");
    double sum=0;
    for(int i=0 ; i<length ; ++i) {
      double LL=fel->logLikelihood(i);
      sum+=LL;
      file.write(sum);
    }
    file.close();
    
    return 0;
  }

