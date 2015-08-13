/****************************************************************
 infer-ancestral-states.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/DnaAlphabet.H"
#include "Phylogeny.H"
#include "FitchParsimony.H"
using namespace std;
using namespace BOOM;


class Application
{
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



Application::Application()
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=3)
      throw String("infer-ancestral-states <in.phy> <in.maf> <out.maf>");
    String phyloFile=cmd.arg(0);
    String mafFile=cmd.arg(1);
    String outFile=cmd.arg(2);

    // Load the phylogeny
    Phylogeny phylogeny(phyloFile);

    // Load the alignments
    MultiAlignment *alignment;
    ifstream is(mafFile.c_str());
    ofstream os(outFile.c_str());
    if(!is.good()) throw String("Can't open file: ")+mafFile;
    Symbol gapSymbol=DnaAlphabet::global.lookup('N');
    while(alignment=MultiAlignment::nextAlignmentFromMAF(is)) {
      alignment->toupper();
      MultSeqAlignment A(*alignment,DnaAlphabet::global,gapSymbol);
      
      // Infer the ancestral states
      FitchParsimony fitch(phylogeny,DnaAlphabet::global,A,gapSymbol);
      fitch.run();
      
      // Append to output file
      A.save(os);
      os<<endl;

      // Clean up 
      delete alignment;
    }
    os.close();

    return 0;
  }

