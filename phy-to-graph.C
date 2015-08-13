/****************************************************************
 phy-to-graph.C
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
    if(cmd.numArgs()!=1)
      throw String("phy-to-graph <*.phy>");
    String infile=cmd.arg(0);
    Phylogeny phy(infile);
    Vector<PhylogenyNode*> nodes;
    phy.gatherNodes(nodes);
    int n=nodes.size();
    for(int i=0 ; i<n ; ++i) {
      PhylogenyNode *node=nodes[i];
      String name=node->getName();
      cout<<node->getNodeType()<<" "<<name<<endl;
      PhylogenyNode *parent=node->getParent();
      if(parent)
	cout<<name<<" -> "<<parent->getName()<<" : "
	    <<node->getDistanceToParent()<<endl;
    }
    
    return 0;
  }

