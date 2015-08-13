/************************************************************************
 UPGMA.C
 Bill Majoros - bmajoros@duke.edu

 Open source; Gnu General Public License (GPL) version 3.
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
  MultiAlignment alignment;
  int nextAncestor;
  Map<PhylogenyNode*,double> cladeHeights;

  PhylogenyNode *buildTree();
  void computePairwiseDistances(TriangMatrix<float> &);
  float computeDistance(int trackI,int trackJ);
  void mergeSubtrees(Vector<PhylogenyNode*> &nodes,
		     TriangMatrix<float> &distanceMatrix,
		     int numSubtrees);
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
  : nextAncestor(1)
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String(
"\nUPGMA <alignment.maf> <outfile>\n\n");
  String infile=cmd.arg(0);
  String outfile=cmd.arg(1);
  
  // Load the alignment
  alignment.loadMAF(infile);

  // Infer the tree
  PhylogenyNode *tree=buildTree();

  // Output the tree
  ofstream os(outfile.c_str());
  tree->output(os);

  return 0;
}



PhylogenyNode *Application::buildTree()
{
  // First, build a matrix of all pairwise distances
  int numTracks=alignment.getNumTracks();
  int alignmentLength=alignment.getLength();
  TriangMatrix<float> distanceMatrix(numTracks);
  computePairwiseDistances(distanceMatrix);

  // Initialize a set of subtrees, one per species
  Vector<PhylogenyNode*> nodes;
  for(int i=0 ; i<numTracks ; ++i)
    {
      LeafNode *node=new LeafNode;
      node->setName(alignment.getIthTrack(i).getName());
      nodes.push_back(node);
      cladeHeights[node]=0;
    }

  // Now iteratively merge subtrees until only one tree remains
  for(int i=numTracks ; i>1 ; --i)
    mergeSubtrees(nodes,distanceMatrix,i);
  return nodes[0];
}



void Application::computePairwiseDistances(
			     TriangMatrix<float> &distanceMatrix)
{
  int numTracks=alignment.getNumTracks();
  int alignmentLength=alignment.getLength();
  for(int i=0 ; i<numTracks ; ++i)
    for(int j=i+1 ; j<numTracks ; ++j)
      distanceMatrix(i,j)=computeDistance(i,j);
}



float Application::computeDistance(int iIndex,int jIndex)
{
  AlignmentTrack &iTrack=alignment.getIthTrack(iIndex);
  AlignmentTrack &jTrack=alignment.getIthTrack(jIndex);
  int iLength=iTrack.getLength();
  int jLength=jTrack.getLength();
  int alignedLength=min(iLength,jLength);
  int matches=0, mismatches=0;
  for(int pos=0 ; pos<alignedLength ; ++pos)
    {
      char iChar=iTrack[pos], jChar=jTrack[pos];
      if(iChar==jChar)
	if(iChar=='-' || iChar==' ') continue;
	else ++matches;
      else ++mismatches;
    }
  float distance=mismatches/float(matches+mismatches);
  return distance;
}



void Application::mergeSubtrees(Vector<PhylogenyNode*> &nodes,
				TriangMatrix<float> &distanceMatrix,
				int numSubtrees)
{
  // Find the closest pair of subtrees
  pair<int,int> bestCell;
  float bestScore=POSITIVE_INFINITY;
  for(int i=0 ; i<numSubtrees ; ++i)
    for(int j=i+1 ; j<numSubtrees ; ++j)
      {
	float score=distanceMatrix(i,j);
	if(score<bestScore) 
	  {
	    bestScore=score;
	    bestCell=pair<int,int>(i,j);
	  }
      }
  
  // Merge these two subtrees by connecting them via a new internal node
  int first=bestCell.first, second=bestCell.second;
  PhylogenyNode *left=nodes[first];
  PhylogenyNode *right=nodes[second];
  InternalNode *newNode=new InternalNode(left,right);
  newNode->setName(String("A")+nextAncestor);
  ++nextAncestor;
  double dist=bestScore-cladeHeights[left]-cladeHeights[right];
  float halfDist=dist/2;
  newNode->setLeftDistance(halfDist);
  newNode->setRightDistance(halfDist);
  cladeHeights[newNode]=halfDist+cladeHeights[left];

  // Replace the first subtree with the newly merged subtree
  for(int i=0 ; i<numSubtrees ; ++i)
    {
      if(i==first || i==second) continue;
      float c1=nodes[first]->getNumLeaves();
      float c2=nodes[second]->getNumLeaves();
      float d1=distanceMatrix.safeIndex(first,i);
      float d2=distanceMatrix.safeIndex(second,i);
      distanceMatrix.safeIndex(first,i)=(c1*d1+c2*d2)/(c1+c2);
    }
  nodes[first]=newNode;

  // Replace the second subtree with the last subtree in the matrix
  int lastIndex=numSubtrees-1;
  if(lastIndex==first || lastIndex==second) return;
  nodes[second]=nodes[lastIndex];
  for(int i=0 ; i<numSubtrees ; ++i)
    {
      if(i==second) continue;
      distanceMatrix.safeIndex(second,i)=
	distanceMatrix.safeIndex(lastIndex,i);
    }
}



