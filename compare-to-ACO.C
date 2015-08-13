/****************************************************************
 compare-to-ACO.C
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
#include "Phylogeny.H"
#include "NthOrdRateMatrix.H"
#include "RateMatrixType.H"
#include "TRCO_Felsenstein.H"
#include "RCO_Felsenstein.H"
#include "ACO_Felsenstein.H"
#include "HOG_Felsenstein.H"
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
    CommandLine cmd(argc,argv,"c:");
    if(cmd.numArgs()!=3)
        throw String(
"\ncompare-to-ACO [-c <context-type>] <*.phy> <*.maf> <*.matrix>\n\
\n\
   where <context-type> = HOG | RCO | TRCO | MP | AMINO | CODON\n\
\n\n");
    String phyFile=cmd.arg(0);
    String mafFile=cmd.arg(1);
    String matrixFile=cmd.arg(2);
    int numThreads=Environment::lookup("NUM_CPUS");
    if(numThreads<1) numThreads=1;
    
    // Load the alignments
    Vector<MultiAlignment*> alignments;
    MultiAlignment::loadMAF(mafFile,alignments);
    MultiAlignment *alignment=MultiAlignment::combine(alignments,true);
    alignment->toupper();
    Symbol gapSymbol=DnaAlphabet::global.lookup('N');
    MultSeqAlignment A(*alignment,DnaAlphabet::global,gapSymbol);

    // Load the phylogeny
    Phylogeny *phylogeny=new Phylogeny(phyFile);

    // Load the rate matrix
    NthOrdRateMatrix *Q=NthOrdRateMatrix::load(matrixFile);
    int order=Q->getOrder();
    bool dual=Q->isDual();

    // Identify the desired context type
    ContextType contextType=CT_MP;
    if(cmd.option('c')) contextType=contextTypeFromString(cmd.optParm('c'));
        
    // Compute likelihood
    double LL=0;
    switch(contextType) {
        case CT_MP: {
            if(dual) {
                Q=Q->ancestralContextsOnly();
                order/=2;
		dual=false;
            }
            FitchParsimony fitch(*phylogeny,DnaAlphabet::global,A,gapSymbol);
            fitch.run();
            phylogeny->attachAlignment(A);
            FitchFelsenstein fel(order,*phylogeny,*Q,A,DNA,numThreads);
            LL=fel.logLikelihood();
        }
            break;
        case CT_TRCO: {
            if(dual) {
                Q=Q->ancestralContextsOnly();
                order/=2;
		dual=false;
            }
            phylogeny->attachAlignment(A);
            TRCO_Felsenstein fel(order,*phylogeny,*Q,A,DNA,numThreads);
            LL=fel.logLikelihood();
        }
            break;
        case CT_RCO: {
            if(dual) {
                Q=Q->ancestralContextsOnly();
                order/=2;
		dual=false;
            }
            phylogeny->attachAlignment(A);
            RCO_Felsenstein fel(order,*phylogeny,*Q,A,DNA,numThreads);
            LL=fel.logLikelihood();
        }
	  break;
        case CT_HOG: {
            if(!dual) throw "HOG model requires dual contexts";
            phylogeny->attachAlignment(A);
            HOG_Felsenstein fel(order/2,*phylogeny,*Q,A,DNA,numThreads);
            LL=fel.logLikelihood();
        }
	  break;
    }

    if(dual) {
      Q=Q->ancestralContextsOnly();
      order/=2;
    }
    phylogeny->attachAlignment(A);
    ACO_Felsenstein acoFel(order,*phylogeny,*Q,A,DNA,numThreads);
    double acoLL=acoFel.logLikelihood();

    // Generate output
    double diff=LL-acoLL;
    cout<<diff<<endl;

    return 0;
  }

