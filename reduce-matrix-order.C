/****************************************************************
 reduce-matrix-order.C
 william.majoros@duke.edu

 This is open-source software,
 governed by the Gnu General Public License (GPL) version 3 (see www.opensource.org).
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
    CommandLine cmd(argc,argv,"ad");
    if(cmd.numArgs()!=3)
      throw String(
"\nreduce-matrix-order [options] <in.matrix> <new-order> <out.matrix>\n\
   options:\n\
      -a = average lower-order models from higher-order models\n\
      -d = use dual contexts (NOTE: order is for parent only!)\n\
\n");
    String infile=cmd.arg(0);
    int newOrder=cmd.arg(1).asInt();
    String outfile=cmd.arg(2);

    // Load matrix
    NthOrdRateMatrix *M=NthOrdRateMatrix::load(infile);
    if(cmd.option('a')) M->averageLowerOrderModel(cmd.option('d'));
    
    // Change order
    if(newOrder<M->getOrder()) {
        if(cmd.option('d')) newOrder*=2;
        NthOrdRateMatrix *newM=M->getLowerOrderModel(newOrder);
        M=newM;
    }

    // Write out new matrix file
    M->save(outfile);
    
    return 0;
  }

