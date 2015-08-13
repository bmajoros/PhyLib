/****************************************************************
 test-felsenstein.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/GSL/Matrix.H"
#include "BOOM/GSL/Vector.H"
#include "BOOM/Random.H"
#include "FelsensteinsAlgorithm.H"
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
    if(cmd.numArgs()!=4)
      throw String("test-felsenstein <alignment.maf> <*.phy> <Q> <rootName>");
    String alignFile=cmd.arg(0);
    String phyFile=cmd.arg(1);
    String matrixFile=cmd.arg(2);
    String rootName=cmd.arg(3);

    // Load data
    Phylogeny phylogeny(phyFile);
    RateMatrix Q(matrixFile,DNA);
    MultiAlignment A;
    A.loadMAF(alignFile);


    //#### perturb
    //Q(4,2)*=2;
    //Q.installDiagonals();
    //####

    // Run the algorithm
    if(!phylogeny.reroot(rootName))
      throw "Error rerooting the tree; possibly the wrong track name?";
    FelsensteinsAlgorithm F(phylogeny,Q,A,DNA);
    double L=F.logLikelihood();
    cout<<"likelihood: "<<L<<endl;
  }
