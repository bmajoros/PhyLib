/****************************************************************
 train-parallel.C
 bmajoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3.
 ****************************************************************/
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/MultiAlignment.H"
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/Exceptions.H"
#include "BOOM/AlphabetMap.H"
#include "BOOM/PureDnaAlphabet.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                       class Application
 ****************************************************************/
class Application
{
  CommandLine *cmd;
  MultSeqAlignment *alignment;
  Alphabet *alphabet;
  BitSet gapSymbols;
  AlphabetMap *alphabetMap; // maps gapped alphabet to ungapped alphabet
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
    catch(const String &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(const RootException &e)
      {
	cerr << "exception: "<< e.getMessage() << endl;
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
    : gapSymbols(256)
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  cmd=new CommandLine(argc,argv,"");
  if(cmd->numArgs()!=2)
    throw String("\ndelete-target-gaps <in.maf> <target-name>\n");
  String mafFile=cmd->arg(0);
  String targetName=cmd->arg(1);
  
  // Setting up alphabet
  alphabet=&DnaDashDotAlphabet::global();
  gapSymbols.addMember(alphabet->lookup('-'));
  gapSymbols.addMember(alphabet->lookup('.'));
  gapSymbols.addMember(alphabet->lookup('N'));
  alphabetMap=new DropGapMapping(*alphabet,PureDnaAlphabet::global());
  
  // Load the alignments
  //int rootID=phylogeny->getRoot()->getID();
  alignment=NULL;
  ifstream is(mafFile.c_str());
  if(!is.good()) throw String("Can't open file ")+mafFile;
  while(!is.eof()) {
    MultiAlignment *a=new MultiAlignment;
    a->loadMAF(is);
    if(a->getNumTracks()>0) {
      a->toupper();
      MultSeqAlignment *m=new MultSeqAlignment(*a,*alphabet,gapSymbols);
      int rootID=m->getTrackByName(targetName).getID();
      m->deleteTargetGaps(rootID);
      if(!alignment) alignment=m;
      else {
	alignment->append(*m);
	delete m;
      }
    }
    delete a;
  }
  is.close();

  // Emit the processed alignment
  alignment->printSlice(cout,0,alignment->getLength(),'+',60);

  return 0;
}
