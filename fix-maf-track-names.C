/************************************************************************
 fix-maf-track-names.C
 Bill Majoros - bmajoros@duke.edu
************************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/MultiAlignment.H"
#include "BOOM/Regex.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/TriangMatrix.H"
#include "BOOM/Constants.H"
#include "BOOM/Map.H"
#include "SubstitutionMatrix.H"
#include "Phylogeny.H"
using namespace std;
using namespace BOOM;


/************************************************************************
 */
class Application
{
public:
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

int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"m:");
  if(cmd.numArgs()!=3)
    throw String(
"\nfix-maf-track-names <in.maf> <new-track-0-name> <out.maf>\n\
\n");
  String infile=cmd.arg(0);
  String targetName=cmd.arg(1);
  String outfile=cmd.arg(2);

  ifstream is(infile.c_str());
  ofstream os(outfile.c_str());
  if(!is.good()) throw infile+" : can't read file";
  while(!is.eof()) {
    MultiAlignment *alignment=MultiAlignment::nextAlignmentFromMAF(is);
    if(!alignment) break;
    alignment->getIthTrack(0).setName(targetName);
    os<<*alignment<<endl;
    delete alignment;
  }

  return 0;
}



