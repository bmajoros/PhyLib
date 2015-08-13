/****************************************************************
 get-phy-distance.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "Phylogeny.H"
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
    throw String("get-phy-distance <*.phy> <leaf1> <leaf2>");
  String infile=cmd.arg(0);
  String name1=cmd.arg(1);
  String name2=cmd.arg(2);
  
  // Load phylogeny
  Phylogeny phy(infile);
  Vector<PhylogenyNode*> *nodes=phy.gatherNodes();

  // Get distance
  PhylogenyNode *node1=phy.findNode(name1), *node2=phy.findNode(name2);
  if(!node1) throw String("Can't find taxon ")+name1;
  if(!node2) throw String("Can't find taxon ")+name2;
  double dist=phy.distanceBetween(node1,node2);
  cout<<dist<<endl;

  return 0;
}

