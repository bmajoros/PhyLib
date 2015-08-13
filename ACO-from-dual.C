/****************************************************************
 ACO-from-dual.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "NthOrdRateMatrix.H"
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
    if(cmd.numArgs()!=2)
      throw String("\nACO-from-dual <in.matrix> <out.matrix>\n");
    String infile=cmd.arg(0);
    String outfile=cmd.arg(1);

    // Load matrix
    NthOrdRateMatrix *M=NthOrdRateMatrix::load(infile);

    // Convert to ACO from dual
    NthOrdRateMatrix *newM=M->ancestralContextsOnly();

    // Write out new matrix file
    newM->save(outfile);
    
    return 0;
  }

