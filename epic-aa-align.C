/****************************************************************
 epic-aa-align.C
 Bill Majoros - bmajoros@duke.edu
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "tigr++/TigrCommandLine.H"
#include "tigr++/TigrFastaReader.H"
#include "tigr++/AminoAlphabet.H"
#include "tigr++/DnaAlphabet.H"
#include "tigr++/MultiAlignment.H"
#include "tigr++/TigrRegex.H"
#include "tigr++/TigrProteinTrans.H"
#include "alignment/MultiHAF.H"
#include "alignment/SubstitutionMatrix.H"
using namespace std;


class Application
{
  TigrRegex lineParser;
  TigrVector<TigrString> labels;

  void emitBlock(MultiHAF<float> &);
  void process(int phase,SubstitutionMatrix<float> &M,float gapPenalty,
	       int bandWidth,const TigrString &infile,
	       float frameshiftPenalty);
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
  : lineParser(
    "s\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(.)\\s+(\\d+)\\s+(\\S+)")
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  TigrCommandLine cmd(argc,argv,"qn:w:r:");
  if(cmd.numArgs()!=5)
    throw TigrString(
"\nepic-aa-align <SubstMatrix> <GapPenalty> <NucAlign.maf> <BandWidth>\n\
              <FrameshiftPenalty>\n\n\
example: epic-aa-align blosum62 5 human.ENm004.maf 5 8\n");
  TigrString matrixFile=cmd.arg(0);
  float gapPenalty=-fabs(cmd.arg(1).asFloat());
  TigrString infile=cmd.arg(2);
  int bandWidth=cmd.arg(3).asInt();
  float frameshiftPenalty=-fabs(cmd.arg(4).asFloat());
  
  // Load the substitution matrix
  Alphabet &alphabet=AminoAlphabet::global;
  SubstitutionMatrix<float> M(matrixFile,alphabet);
  
  // Process the nuc multi-alignment file
  for(int phase=0 ; phase<3 ; ++phase)
    process(phase,M,gapPenalty,bandWidth,infile,frameshiftPenalty);

  return 0;
}



void Application::process(int phase,
			  SubstitutionMatrix<float> &M,
			  float gapPenalty,
			  int bandWidth,
			  const TigrString &infile,
			  float frameshiftPenalty)
{
  TigrVector<Sequence*> seqs;
  MultiHAF<float> aligner(M,gapPenalty,bandWidth,frameshiftPenalty);
  ifstream is(infile.c_str());
  if(!is.good()) throw TigrString("Can't open ")+infile;
  while(!is.eof())
    {
      TigrString line;
      line.getline(is);
      if(is.eof() && line.isEmpty()) break;
      line.trimWhitespace();
      if(line.isEmpty()) 
	{
	  if(seqs.size()==0) continue;
	  aligner.performAlignment();
	  emitBlock(aligner);
	  aligner.reset();
	  int n=seqs.size();
	  for(int i=0 ; i<n ; ++i) delete seqs[i];
	  seqs.clear();
	  labels.clear();
	  continue;
	}
      if(line[0]=='s')
	{
	  if(lineParser.match(line))
	    {
	      TigrString name=lineParser[1];
	      int begin=lineParser[2].asInt();
	      int len=lineParser[3].asInt();
	      char strand=lineParser[4][0];
	      int contigLength=lineParser[5].asInt();
	      TigrString seqStr=lineParser[6];
	      seqStr.toupper();
	      seqStr=seqStr.substitute("-","");
	      if(seqs.size()==0 && phase!=0)
		if(signed(seqStr.length())-phase<3) 
		  continue;
		else
		  seqStr=seqStr.substring(phase,seqStr.length()-phase);
	      Sequence *seq=new Sequence(seqStr,DnaAlphabet::global);
	      seqs.push_back(seq);
	      aligner.addSequence(*seq);
	      labels.push_back(name);
	    }
	  else throw TigrString("Can't parse line: ")+line;
	}
    }
  cout<<"# END OF PHASE "<<phase<<endl;
}



void Application::emitBlock(MultiHAF<float> &aligner)
{
  int numTracks=aligner.numTracks(), longestLabel=0;
  for(int i=0 ; i<numTracks ; ++i)
    {
      int length=labels[i].length();
      if(length>longestLabel) longestLabel=length;
    }
  for(int i=0 ; i<numTracks ; ++i)
    {
      TigrString label=labels[i];
      TigrString padding(longestLabel-label.length(),' ');
      label+=padding;
      const TigrString &track=aligner.getIthTrack(i);
      cout<<"s "<<label<<" 0 0 . 0 "<<track<<endl;
    }
  cout<<endl;
}
