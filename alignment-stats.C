/************************************************************************
 alignment-stats.C
 bmajoros@duke.edu

 Open source software under Gnu General Public License (GPL) version 3.
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

struct TrackInfo {
    TrackInfo() : numGaps(0), numGC(0), matchTarget(0), noGaps(0) {}
    int numGaps;
    int noGaps; // #cols in which this track & the target are both non-gap
    int numGC;
    int matchTarget;
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
  if(cmd.numArgs()!=1)
    throw String("\nalignment-stats <alignment.maf>\n\n");
  String infile=cmd.arg(0);
  
  // Process the alignments
  ifstream is(infile.c_str());
  if(!is.good()) throw infile+" : can't open file";
  int totalLength=0, totalDNA=0, numTracks=0, numAlignments=0;
  Map<String,TrackInfo> trackInfo;
  while(true) {
      MultiAlignment *alignment=MultiAlignment::nextAlignmentFromMAF(is);
      if(!alignment) break;
      ++numAlignments;
      alignment->toupper();
      totalLength+=alignment->getLength();
      int n=alignment->getNumTracks();
      totalDNA+=alignment->getLength()*n;
      if(n>numTracks) numTracks=n;
      AlignmentTrack &targetTrack=alignment->getIthTrack(0);
      for(int i=0 ; i<n ; ++i) {
          AlignmentTrack &track=alignment->getIthTrack(i);
          String name=track.getName();
          if(i==0) name="(target)";
          int L=track.getLength();
          TrackInfo &info=trackInfo[name];
          for(int j=0 ; j<L ; ++j) {
              char c=track[j], t=targetTrack[j];
              if(c=='-' || c=='.') ++info.numGaps;
              if(c=='G' || c=='C') ++info.numGC;
              if(c==targetTrack[j] &&
                 (c=='A'||c=='C'||c=='G'||c=='T')) ++info.matchTarget;
              if((c=='A'||c=='C'||c=='G'||c=='T') &&
                 (t=='A'||t=='C'||t=='G'||t=='T')) ++info.noGaps;
          }

      }
      delete alignment;
  }
  cout<<"Total length: "<<totalLength<<", total DNA: "<<totalDNA<<endl;
  cout<<"Number of tracks: "<<numTracks<<endl;
  cout<<"Number of alignments: "<<numAlignments<<endl;
  
  Map<String,TrackInfo>::iterator cur=trackInfo.begin(), end=trackInfo.end();
  for(; cur!=end ; ++cur) {
      String name=(*cur).first;
      TrackInfo &track=(*cur).second;
      int L=totalLength, numGaps=track.numGaps, numGC=track.numGC;
      int numACGT=L-numGaps;
      float gc=int((numGC/float(numACGT))*100+5/9.0);
      float coverage=int(numACGT/float(L)*100+5/9.0);
      float matchTarget=int(track.matchTarget/float(track.noGaps)*100+5/9.0);
      cout<<"Track: "<<name<<" GC="<<gc<<"% coverage="<<coverage<<"% match="
          <<matchTarget<<"%"<<endl;
  }
  
  return 0;
}



