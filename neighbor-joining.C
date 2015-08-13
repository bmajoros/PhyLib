/************************************************************************
 neighbor-joining.C
 see Durbin,Krogh,Mitchison,Eddy, p171

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
  int nextAncestor, numTracks;
  TriangMatrix<float> distanceMatrix;
  TriangMatrix<long int> matchesMatrix, mismatchesMatrix;
  Vector<PhylogenyNode*> nodes;
  Vector<String> trackNames;
  Map<String,int> nameToID;

  PhylogenyNode *buildTree(String infile);
  void computePairwiseDistances(String infile);
  void computeDistance(int trackI,int trackJ,MultiAlignment *);
  void mergeSubtrees(int numSubtrees);
  float ri(int i,int j);
  float Dij(int i,int j);
  void countTracks(String infile);
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
  : nextAncestor(1), numTracks(0)
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"m:");
  if(cmd.numArgs()!=2)
    throw String(
"\nneighbor-joining [-m] <alignment.maf> <outfile>\n\
where: -m <filename> = write distance matrix into file\n\
\n");
  String infile=cmd.arg(0);
  String outfile=cmd.arg(1);
  
  // Infer the tree
  PhylogenyNode *tree=buildTree(infile);

  // Output the tree
  ofstream os(outfile.c_str());
  tree->output(os);

  if(cmd.option('m')) {
    ofstream os(cmd.optParm('m').c_str());
    if(!os.good()) throw cmd.optParm('m')+" : can't write to file";
    os<<distanceMatrix<<endl;
    cout<<matchesMatrix<<"\n\n"<<mismatchesMatrix<<endl;
  }

  return 0;
}



PhylogenyNode *Application::buildTree(String infile)
{
  // First, build a matrix of all pairwise distances
  computePairwiseDistances(infile);

  // Initialize a set of subtrees, one per species
  for(int i=0 ; i<numTracks ; ++i)
    {
      LeafNode *node=new LeafNode;
      node->setName(trackNames[i]);
      nodes.push_back(node);
    }

  // Now iteratively merge subtrees until only two trees remain
  for(int i=numTracks ; i>2 ; --i)
    mergeSubtrees(i);

  // Merge the final two subtrees
  InternalNode *root=new InternalNode(nodes[0],nodes[1]);
  root->setName("root");
  double halfDist=distanceMatrix(0,1)/2;
  root->setLeftDistance(halfDist);
  root->setRightDistance(halfDist);
  return root;
}



void Application::countTracks(String infile) {
  ifstream is(infile.c_str());
  if(!is.good()) throw infile+" : can't read file";
  while(!is.eof()) {
    MultiAlignment *alignment=MultiAlignment::nextAlignmentFromMAF(is);
    if(!alignment) break;
    int n=alignment->getNumTracks();
    for(int i=0 ; i<n ; ++i) {
      AlignmentTrack &track=alignment->getIthTrack(i);
      String name=track.getName();
      if(!nameToID.isDefined(name)) {
	nameToID[name]=numTracks++;
	trackNames.push_back(name);
      }
    }
    delete alignment;
  }
}



void Application::computePairwiseDistances(String infile)
{
  // Initialize sizes of data structures
  countTracks(infile);
  distanceMatrix.resize(numTracks);
  matchesMatrix.resize(numTracks);
  mismatchesMatrix.resize(numTracks);
  distanceMatrix.setAllTo(0);
  matchesMatrix.setAllTo(0);
  mismatchesMatrix.setAllTo(0);

  // Load the alignments
  ifstream is(infile.c_str());
  if(!is.good()) throw infile+" : can't read file";
  while(!is.eof()) {
    MultiAlignment *alignment=MultiAlignment::nextAlignmentFromMAF(is);
    if(!alignment) break;
    alignment->toupper();
    for(int i=0 ; i<numTracks ; ++i)
      for(int j=i+1 ; j<numTracks ; ++j) {
	computeDistance(i,j,alignment);
      }
    delete alignment;
  }
  for(int i=0 ; i<numTracks ; ++i)
    for(int j=i+1 ; j<numTracks ; ++j) {
      int matches=matchesMatrix(i,j), mismatches=mismatchesMatrix(i,j);
      float distance=mismatches/float(matches+mismatches);
      distanceMatrix(i,j)=distance;
    }
}



void Application::computeDistance(int iIndex,int jIndex,
				  MultiAlignment *alignment)
{
  AlignmentTrack &iTrack=alignment->getTrackByName(trackNames[iIndex]);
  AlignmentTrack &jTrack=alignment->getTrackByName(trackNames[jIndex]);
  int iLength=iTrack.getLength();
  int jLength=jTrack.getLength();
  int alignedLength=min(iLength,jLength);
  int matches=0, mismatches=0;
  for(int pos=0 ; pos<alignedLength ; ++pos)
    {
      char iChar=iTrack[pos], jChar=jTrack[pos];
      if(iChar=='-' || iChar==' ' || jChar=='-' || jChar==' ') continue;
      if(iChar==jChar) ++matches;
      else ++mismatches;
    }
  matchesMatrix(iIndex,jIndex)+=matches;
  mismatchesMatrix(iIndex,jIndex)+=mismatches;
}



void Application::mergeSubtrees(int numSubtrees)
{
  // Find the closest pair of subtrees
  pair<int,int> bestCell;
  float bestScore=POSITIVE_INFINITY;
  for(int i=0 ; i<numSubtrees ; ++i)
    for(int j=i+1 ; j<numSubtrees ; ++j)
      {
	float score=Dij(i,j);
	if(score<bestScore) 
	  {
	    bestScore=score;
	    bestCell=pair<int,int>(i,j);
	  }
      }
  
  // Merge these two subtrees by connecting them via a new internal node
  int first=bestCell.first, second=bestCell.second;
  const double distIJ=distanceMatrix(first,second);
  PhylogenyNode *left=nodes[first];
  PhylogenyNode *right=nodes[second];
  InternalNode *newNode=new InternalNode(left,right);
  newNode->setName(String("A")+nextAncestor);
  ++nextAncestor;

  double leftDistance=(distIJ+ri(first,second)-ri(second,first))/2;
  double rightDistance=distIJ-leftDistance;
  newNode->setLeftDistance(leftDistance);
  newNode->setRightDistance(rightDistance);

  // Replace the first subtree with the newly merged subtree
  for(int i=0 ; i<numSubtrees ; ++i)
    {
      if(i==first || i==second) continue;
      float d1=distanceMatrix.safeIndex(first,i);
      float d2=distanceMatrix.safeIndex(second,i);
      distanceMatrix.safeIndex(first,i)=(d1+d2-distIJ)/2;
    }
  nodes[first]=newNode;

  // Replace the second subtree with the last subtree in the matrix,
  // and copy its distances over as well (since its index is changing)
  int lastIndex=numSubtrees-1;
  if(lastIndex!=second)
    {
      nodes[second]=nodes[lastIndex];
      for(int i=0 ; i<numSubtrees-1 ; ++i)
	{
	  if(i==second) continue;
	  distanceMatrix.safeIndex(second,i)=
	    distanceMatrix.safeIndex(lastIndex,i);
	}
    }
  nodes.resize(numSubtrees-1);
}



float Application::ri(int i,int j)
{
  double sum=0;
  int n=nodes.size();
  for(int k=0 ; k<n ; ++k)
    if(k!=i)
      sum+=distanceMatrix.safeIndex(i,k);
  double ave=sum/(n-2);
  return ave;
}



float Application::Dij(int i,int j)
{
  return distanceMatrix.safeIndex(i,j)-ri(i,j)-ri(j,i);
}


