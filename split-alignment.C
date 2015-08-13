/************************************************************************
 split-alignment.C
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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String(
"\nsplit-alignment <alignment.maf> <output-directory>\n\n");
  String infile=cmd.arg(0);
  String dir=cmd.arg(1);

  // Process the alignments
  int fileNum=1;
  ifstream is(infile.c_str());
  if(!is.good()) throw infile+" : can't open file";
  while(true) {
      MultiAlignment *alignment=MultiAlignment::nextAlignmentFromMAF(is);
      if(!alignment) break;
      String outfile=dir+"/"+fileNum+".maf";
      cout<<outfile<<endl;
      ofstream os(outfile.c_str());
      if(!os.good()) throw outfile+" : can't create file";
      os<<*alignment<<endl;
      os.close();
      delete alignment;
      ++fileNum;
  }
  
  return 0;
}



