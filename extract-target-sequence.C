/****************************************************************
 extract-target-sequence.C
 Bill Majoros - bmajoros@duke.edu 2/9/05

 This is open-source software  governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/MultiAlignment.H"
#include "BOOM/FastaWriter.H"
using namespace std;
using namespace BOOM;

class Application
{
public:
  int main(int argc,char *argv[]);
};


/**************************************************************************
 main() -- just handles exceptions
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



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw string("extract-target-sequence <*.maf> <ID> <outfile.fasta>");
  String alignmentFile=cmd.arg(0);
  int ID=cmd.arg(1).asInt();
  String outfile=cmd.arg(2);
  
  // Process each alignment in the input file
  FastaWriter writer;
  ifstream is(alignmentFile.c_str());
  ofstream os(outfile.c_str());
  if(!is.good()) throw String("Error opening file ")+alignmentFile;
  while(!is.eof()) {
    MultiAlignment *alignment=MultiAlignment::nextAlignmentFromMAF(is);

    // Extract zeroth track, remove dashes, convert to uppercase
    AlignmentTrack &targetTrack=alignment->getIthTrack(0);
    String seqStr=targetTrack.getSeq().substitute("-","");
    seqStr.toupper();
    
    // Write into FASTA file
    String defline=String(">")+ID;
    writer.addToFasta(defline,seqStr,os);

    delete alignment;
    ++ID;
  }
  is.close();
  os.close();
    
  return 0;
}




