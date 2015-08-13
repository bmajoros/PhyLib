/****************************************************************
 make-phase-0.C
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
#include "BOOM/Vector.H"
#include "BOOM/Random.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/Exceptions.H"
#include "BOOM/BitSet.H"
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
  void processAlignments(const String &infile,const String &outfile);
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



Application::Application()
    : alignment(NULL), gapSymbols(256)
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  cmd=new CommandLine(argc,argv,"");
  if(cmd->numArgs()!=2)
    throw String("\nmake-phase-0 <in.maf> <out.maf>\n\n");
  String infile=cmd->arg(0);
  String outfile=cmd->arg(1);

  // Setting up alphabet
  alphabet=&DnaDashDotAlphabet::global();
  gapSymbols.addMember(alphabet->lookup('-'));
  gapSymbols.addMember(alphabet->lookup('.'));
  gapSymbols.addMember(alphabet->lookup('N'));

  // Load the alignments
  processAlignments(infile,outfile);
  return 0;
}



void Application::processAlignments(const String &infile,const String &outfile)
{
  alignment=NULL;
  ifstream is(infile.c_str());
  if(!is.good()) throw String("Can't open file ")+infile;
  while(!is.eof()) {
    MultiAlignment *a=new MultiAlignment;
    a->loadMAF(is);
    if(a->getNumTracks()>0 && a->getLength()>0) {
      a->toupper();
      MultSeqAlignment *m=new MultSeqAlignment(*a,*alphabet,gapSymbols);
      //m->deleteTargetGaps(phylogeny->getRoot()->getID());
      int phase=a->getPhase();
      int len=a->getLength();
      int begin=(3-phase)%3;
      int end=len;
      int endPhase=(phase+len)%3;
      end-=endPhase;
      int sliceLen=end-begin;
      if(sliceLen%3>0) throw "Application::loadAlignments";//###DEBUGGING
      MultSeqAlignment *slice=m->getSlice(begin,sliceLen);
      delete m;
      if(!alignment) alignment=slice;
      else {
	alignment->append(*slice);
	delete slice;
      }
    }
    delete a;
  }
  is.close();
  alignment->setPhase(0);
  alignment->save(outfile);
}
