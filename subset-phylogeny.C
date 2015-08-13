/****************************************************************
 subset-phylogeny.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Set.H"
#include "Phylogeny.H"
using namespace std;
using namespace BOOM;


class Application
{
  Set<String> keep; // taxa
  void parseSpeciesList(const String &speciesList,Set<String> &into);
  void prunePhylogeny(Phylogeny &);
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
    throw String("subset-phylogeny <in.phy> species1,species2,... <out.phy>");
  String infile=cmd.arg(0), speciesList=cmd.arg(1), outfile=cmd.arg(2);

  // Load the old phylogeny
  Phylogeny oldPhy(infile);

  // Prune away unwanted taxa
  parseSpeciesList(speciesList,keep);
  if(keep.size()<2) throw "specify at least two species";
  prunePhylogeny(oldPhy);

  // Save the new phylogeny
  oldPhy.save(outfile);

  return 0;
}



void Application::parseSpeciesList(const String &speciesList,
				   Set<String> &S)
{
  Vector<String> *fields=speciesList.getFields(",");
  Vector<String>::iterator cur=fields->begin(), end=fields->end();
  for(; cur!=end ; ++cur) S.insert(*cur);
  delete fields;
}



void Application::prunePhylogeny(Phylogeny &phy)
{
  phy.prune(keep);
}


