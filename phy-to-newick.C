/****************************************************************
 phy-to-newick.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "../PhyLib/Phylogeny.H"
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
      throw String("phy-to-newick <infile> <outfile>");
    String infile=cmd.arg(0);
    String outfile=cmd.arg(1);

    Phylogeny phylogeny(infile);
    ofstream os(outfile.c_str());
    os<<"#NEXUS\n\nBegin trees;\n\ntree Phylogeny = ";
    phylogeny.printNewick(os);
    os<<";\n\nEnd;\n\n";
    os.close();

    return 0;
  }

