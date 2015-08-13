/****************************************************************
 compare-alignments.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/MultiAlignment.H"
#include "BOOM/Vector.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/MultSeqAlignment.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                       class Application
 ****************************************************************/
class Application
{
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




/****************************************************************
                      Application methods
 ****************************************************************/



Application::Application()
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=2)
        throw String("\ncompare-alignments <1.maf> <2.maf>\n");
    String maf1=cmd.arg(0);
    String maf2=cmd.arg(1);
    
    // Load the alignments
    Vector<MultiAlignment*> alignments1;
    MultiAlignment::loadMAF(maf1,alignments1);
    MultiAlignment *alignment1=MultiAlignment::combine(alignments1,true);
    alignment1->toupper();
    MultSeqAlignment A1(*alignment1,DnaAlphabet::global,
                       DnaAlphabet::global.lookup('N'));

    Vector<MultiAlignment*> alignments2;
    MultiAlignment::loadMAF(maf2,alignments2);
    MultiAlignment *alignment2=MultiAlignment::combine(alignments2,true);
    alignment2->toupper();
    MultSeqAlignment A2(*alignment2,DnaAlphabet::global,
                       DnaAlphabet::global.lookup('N'));

    // Iterate over tracks in one alignment
    int n=A1.getNumTracks();
    for(int i=0 ; i<n ; ++i) {
        AlignmentSeq &track1=A1.getIthTrack(i);
        AlignmentSeq &track2=A2.getTrackByName(track1.getName());
        int L=track1.getLength();
        int correct=0;
        for(int j=0 ; j<L ; ++j)
            if(track1[j]==track2[j]) ++correct;
        float score=correct/float(L);
        score=int(1000*score+5/9.0)/10.0;
        cout<<track1.getName()<<" : "<<score<<"%"<<endl;
    }
    

    return 0;
  }

