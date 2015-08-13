/****************************************************************
 filter-alignment-columns.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/BitSet.H"
#include "BOOM/MultiAlignment.H"
using namespace std;
using namespace BOOM;


class Application
{
  int numTracks;
  bool speciesPredicate(MultiAlignment *,int col,Vector<String> *species);
  bool numPredicate(MultiAlignment *,int col,int minCount);
  MultiAlignment *filter(MultiAlignment *oldAlign,const BitSet &goodCols);
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
  CommandLine cmd(argc,argv,"n:s:");
  if(cmd.numArgs()!=2)
    throw String("\nfilter-alignment-columns [options] <in.maf> <out.maf>\n\
where options are:\n\
    -n # : discard columns having fewer than # non-gap characters\n\
    -s x,y,z : discard any column having a gap in the x, y, or z track\n\
\n");
  String infile=cmd.arg(0), outfile=cmd.arg(1);
  bool filterByNumber=cmd.option('n');
  bool filterBySpecies=cmd.option('s');
  int minCount=0;
  Vector<String> *requiredSpecies=NULL;
  if(filterByNumber) minCount=cmd.optParm('n').asInt();
  if(filterBySpecies) requiredSpecies=cmd.optParm('s').getFields(",");
  
  ifstream is(infile.c_str());
  if(!is.good()) throw String("error opening file ")+infile;
  ofstream os(outfile.c_str());
  if(!os.good()) throw String("error writing to file ")+outfile;
  while(!is.eof()) {
    MultiAlignment *alignment=MultiAlignment::nextAlignmentFromMAF(is);
    if(!alignment) break;
    int numCols=alignment->getLength();
    numTracks=alignment->getNumTracks();
    BitSet goodColumns(numCols);
    for(int col=0 ; col<numCols ; ++col) {
      bool keepColumn=true;
      if(filterByNumber) 
	keepColumn=keepColumn && numPredicate(alignment,col,minCount);
      if(filterBySpecies) 
	keepColumn=keepColumn && speciesPredicate(alignment,col,
						  requiredSpecies);
      if(keepColumn) goodColumns.addMember(col);
    }
    MultiAlignment *filteredAlignment=filter(alignment,goodColumns);
    os<<*filteredAlignment<<endl;
    delete alignment;
    delete filteredAlignment;
  }
  return 0;
}



bool Application::speciesPredicate(MultiAlignment *alignment,int col,
				   Vector<String> *requiredSpecies)
{
  Vector<String>::iterator cur=requiredSpecies->begin(),
    end=requiredSpecies->end();
  for(; cur!=end ; ++cur) {
    const String &species=*cur;
    AlignmentTrack &track=alignment->findOrCreateTrack(species);
    if(track.getLength()<col+1) return false;
    char c=track[col];
    if(c=='.' || c=='-' || c==' ') return false;
  }
  return true;
}



bool Application::numPredicate(MultiAlignment *alignment,int col,int minCount)
{
  int count=0;
  for(int i=0 ; i<numTracks ; ++i) {
    AlignmentTrack &track=alignment->getIthTrack(i);
    if(track.getLength()>col) {
      char c=track[col];
      if(c!='.' && c!='-' && c!=' ') ++count;
    }
  }
  return count>=minCount;
}



MultiAlignment *Application::filter(MultiAlignment *oldAlignment,
				    const BitSet &goodColumns)
{
  int oldLength=oldAlignment->getLength();
  int newLength=goodColumns.cardinality();
  MultiAlignment *newAlignment=new MultiAlignment();
  for(int i=0 ; i<numTracks ; ++i) {
    AlignmentTrack &oldTrack=oldAlignment->getIthTrack(i);
    const String &name=oldTrack.getName();
    AlignmentTrack &newTrack=newAlignment->findOrCreateTrack(name);
    newTrack.extendToLength(newLength);
    int nextPos=0;
    for(int j=0 ; j<oldLength ; ++j) 
      if(goodColumns.isMember(j)) {
	newTrack[nextPos]=oldTrack[j];
	++nextPos;
      }
  }
  return newAlignment;
}

