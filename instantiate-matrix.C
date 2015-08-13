/****************************************************************
 instantiate-matrix.C
 bmajoros@duke.edu
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/DnaAlphabet.H"
#include "NthOrdRateMatrix.H"
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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("\ninstantiate-matrix <infile> <branch-length>\n");
  String infile=cmd.arg(0);
  double branchLength=cmd.arg(1).asFloat();
  
  NthOrdRateMatrix *Q=NthOrdRateMatrix::load(infile);
  NthOrdSubstMatrix *Pt=Q->instantiate(branchLength);
  cout<<*Pt<<endl;

  return 0;
}



