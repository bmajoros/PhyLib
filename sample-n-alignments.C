/************************************************************************
 sample-n-alignments.C
 bmajoros@duke.edu
************************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/MultiAlignment.H"
#include "BOOM/Regex.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/Constants.H"
#include "BOOM/Map.H"
#include "BOOM/File.H"
#include "BOOM/BitSet.H"
#include "BOOM/Random.H"
using namespace std;
using namespace BOOM;


/************************************************************************
 */
class Application
{
  MultiAlignment *alignment;
  int countAlignments(String filename);
  BitSet *getSampleVector(int numHave,int numWant);
public:
  Application();
  int main(int argc,char *argv[]);
};


/************************************************************************
 */
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



/************************************************************************
 */
Application::Application()
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"t:r:");
  if(cmd.numArgs()!=3)
    throw String(
"\nsample-n-alignments [options] <alignment.maf> <n> <out.maf>\n\
   where: -t <trackname> = rename target tracks\n\
          -r <remfile> = put remaining alignments in <remfile>\n\
\n");
  String infile=cmd.arg(0);
  int N=cmd.arg(1);
  String outfile1=cmd.arg(2);
  String outfile2;
  bool wantRemainder=cmd.option('r');
  if(wantRemainder) outfile2=cmd.optParm('r');
  randomize();

  // Create output files
  ofstream os1(outfile1.c_str()), os2;
  if(wantRemainder) os2.open(outfile2.c_str());
  
  // Prepare for sampling
  int n=countAlignments(infile);
  BitSet *sample=getSampleVector(n,N);

  // Process the alignments
  ifstream is(infile.c_str());
  if(!is.good()) throw infile+" : can't open file";
  int index=0;
  while(true) {
    MultiAlignment *alignment=MultiAlignment::nextAlignmentFromMAF(is);
    if(!alignment) break;
    if(cmd.option('t'))
      alignment->getIthTrack(0).rename(cmd.optParm('t'));
    if(sample->isMember(index)) {alignment->printOn(os1);os1<<endl;}
    else if(wantRemainder) {alignment->printOn(os2);os2<<endl;}
    ++index;
  }
  
  return 0;
}



int Application::countAlignments(String filename)
{
  int count=0;
  File file(filename);
  while(!file.eof()) {
    String line=file.getline();
    if(line.substring(0,2)=="a ") ++count;
  }
  return count;
}



BitSet *Application::getSampleVector(int numHave,int numWant) 
{
  // Set first N bits 1s, leave the rest 0s
  BitSet &bs=*new BitSet(numHave);
  for(int i=0 ; i<numWant ; ++i) bs.addMember(i);

  // Now randomly permute the list
  for(int i=0 ; i<numHave ; ++i) 
    bs.swapBits(i,RandomNumber(numHave));

  // Return the list
  return &bs;
}




