/************************************************************************
 first-n-alignments.C
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
using namespace std;
using namespace BOOM;


/************************************************************************
 */
class Application
{
  MultiAlignment *alignment;
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
"\nfirst-n-alignments [options] <alignment.maf> <n> <out.maf>\n\
   where: -t <trackname> = rename target tracks\n\
          -r <remfile> = put remaining alignments in <remfile>\n\
\n");
  String infile=cmd.arg(0);
  int N=cmd.arg(1);
  String outfile1=cmd.arg(2);
  String outfile2;
  bool wantRemainder=cmd.option('r');
  if(wantRemainder) outfile2=cmd.optParm('r');

  // Create output files
  ofstream os1(outfile1.c_str()), os2;
  if(wantRemainder) os2.open(outfile2.c_str());
  ofstream *os=&os1;
  
  // Process the alignments
  ifstream is(infile.c_str());
  if(!is.good()) throw infile+" : can't open file";
  int count=0;
  while(true) {
      MultiAlignment *alignment=MultiAlignment::nextAlignmentFromMAF(is);
      if(!alignment) break;
      if(cmd.option('t'))
          alignment->getIthTrack(0).rename(cmd.optParm('t'));
      if(count<=N || wantRemainder) alignment->printOn(*os);
      delete alignment;
      ++count;
      if(count==N) {
          os=&os2;
          if(!wantRemainder) break;
      }
  }
  int count1=count>N ? N : count;
  int count2=wantRemainder ? count-count1 : 0;
  cout<<count1<<" in file #1, "<<count2<<" in file #2"<<endl;
  
  return 0;
}



