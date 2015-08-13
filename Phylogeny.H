/****************************************************************
 Phylogeny.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Phylogeny_H
#define INCL_Phylogeny_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Vector.H"
#include "BOOM/Array2D.H"
#include "BOOM/Exceptions.H"
#include "ResidueCounter.H"
#include "NthOrdSubstMatrix.H"
#include "IndelHistory.H"
using namespace std;
using namespace BOOM;

class TreeVisitor;

enum BranchEnd {
  PARENT,
  CHILD
};
inline BranchEnd operator!(BranchEnd e) {return e==PARENT ? CHILD : PARENT;}

enum PhylogenyNodeType {
  ROOT_NODE,
  LEAF_NODE,
  INTERNAL_NODE
};

ostream &operator<<(ostream &,PhylogenyNodeType);

enum TraversalOrder {
  PREORDER,
  POSTORDER,
  INORDER // ### not currently supported
};

enum WhichChild {UNIQUE=0,LEFT=0,RIGHT=1};
inline WhichChild otherChild(WhichChild);


enum InsideOutside {
  INSIDE=0,
  OUTSIDE=1
};


/************************************************************************
                        class PhylogenyDecoration
*************************************************************************/
struct PhylogenyDecoration {
  PhylogenyDecoration *clone() {return new PhylogenyDecoration;}
};


/************************************************************************
                          class BranchDecoration
*************************************************************************/
struct BranchDecoration {
  BranchDecoration *clone() {return new BranchDecoration;}
};


/************************************************************************
                            class PhylogenyNode
*************************************************************************/
class PhylogenyBranch;
class PhylogenyNode
{
protected:
  int numLeaves;
  PhylogenyNodeType nodeType;
  PhylogenyNode *parent;
  ResidueCounter residueCounter;
  String name;
  int ID;
  IndelHistory indelHistory;
  PhylogenyDecoration *decoration;
  BitSet cladeMembers;
  Array1D<PhylogenyBranch> branches;
  void setNumLeaves(int);
  PhylogenyNode(const PhylogenyNode &);
public:
  PhylogenyNode(int numLeaves,PhylogenyNodeType);
  virtual ~PhylogenyNode() {}
  BitSet &getCladeMembers() {return cladeMembers;}
  void setName(const String &);
  String &getName();
  void setID(int);
  int getID() const;
  inline PhylogenyNodeType getNodeType() const;
  virtual void output(ostream &)=0;
  virtual void printNewick(ostream &)=0;
  int getNumLeaves() const;
  static PhylogenyNode *load(istream &);
  virtual void postorderTraversal(TreeVisitor &)=0;
  virtual void preorderTraversal(TreeVisitor &)=0;
  virtual void recomputeNumLeaves()=0;
  virtual bool reroot(const String &rootName) {}
  inline PhylogenyNode *getParent() const;
  float getDistanceToParent() const;
  void setParent(PhylogenyNode *);
  inline ResidueCounter &getCounter();
  virtual PhylogenyNode *clone()=0;
  inline IndelHistory &getIndelHistory();
  inline PhylogenyDecoration *&getDecoration();
  WhichChild whichChild();
  void allocateBranches(int);
  PhylogenyBranch &getBranch(int);
  virtual PhylogenyNode *getChild(int i)=0;
  void getPathToRoot(Vector<PhylogenyNode*> &);
};


/************************************************************************
                         class PhylogenyBranch
*************************************************************************/
class PhylogenyBranch {
public:
  PhylogenyBranch(PhylogenyNode *parent,PhylogenyNode *child,int branchID,
		  double length);
  PhylogenyBranch();
  PhylogenyNode *getParent() const;
  PhylogenyNode *getChild() const;
  int getBranchID() const;
  void setBranchID(int id) {branchID=id;}
  double getLength() const;
  void changeLength(double newLength);
  BranchDecoration *&getDecoration();
  void setParent(PhylogenyNode *node) {parent=node;}
  void setChild(PhylogenyNode *node) {child=node;}
private:
  PhylogenyNode *parent, *child;
  int branchID;
  double length;
  BranchDecoration *decoration;
};


/************************************************************************
                            class LeafNode
*************************************************************************/
class LeafNode : public PhylogenyNode {
public:
  LeafNode();
  LeafNode(istream &);
  void output(ostream &);
  virtual void printNewick(ostream &);
  virtual void postorderTraversal(TreeVisitor &);
  virtual void preorderTraversal(TreeVisitor &);
  virtual void recomputeNumLeaves();
  virtual bool reroot(const String &rootName);
  virtual PhylogenyNode *clone();
  virtual PhylogenyNode *getChild(int i) {INTERNAL_ERROR;}
};



/************************************************************************
                            class RootNode
*************************************************************************/
class RootNode : public PhylogenyNode {
  float branchLength;
  PhylogenyNode *child;
  NthOrdSubstMatrix *Pt;
public:
  RootNode(LeafNode &,float branchLength);
  RootNode(float branchLength);
  RootNode(istream &);
  ~RootNode(); // deletes Pt
  void setSubstMatrix(NthOrdSubstMatrix *); // gets deleted!
  NthOrdSubstMatrix *getSubstMatrix();
  void output(ostream &);
  virtual void printNewick(ostream &);
  virtual void postorderTraversal(TreeVisitor &);
  virtual void preorderTraversal(TreeVisitor &);
  PhylogenyNode *getChild();
  void setChild(PhylogenyNode *);
  virtual void recomputeNumLeaves();
  virtual PhylogenyNode *clone();
  float getBranchLength() const;
  void setBranchLength(float);
  virtual PhylogenyNode *getChild(int i) {return child;}
};



/************************************************************************
                           class InternalNode
*************************************************************************/
class InternalNode : public PhylogenyNode {
  PhylogenyNode *left, *right;
  float leftDistance, rightDistance;
  NthOrdSubstMatrix *leftMatrix, *rightMatrix;
  bool resolvedBit; // used by the MUD algorithm
  
  void reroot(PhylogenyNode*,WhichChild,float branchLength);
  PhylogenyNode *getOtherChild(WhichChild);
  float &getDistance(WhichChild);
public:
  InternalNode(PhylogenyNode *left,PhylogenyNode *right);
  InternalNode(istream &);
  ~InternalNode();
  void setLeftSubstMatrix(NthOrdSubstMatrix *); // gets deleted!
  void setRightSubstMatrix(NthOrdSubstMatrix *); // gets deleted!
  NthOrdSubstMatrix *getLeftSubstMatrix();
  NthOrdSubstMatrix *getRightSubstMatrix();
  void output(ostream &);
  virtual void printNewick(ostream &);
  inline PhylogenyNode *getLeft();
  inline PhylogenyNode *getRight();
  int numChildren();
  inline void setLeft(PhylogenyNode *);
  inline void setRight(PhylogenyNode *);
  inline float getLeftDistance() const;
  inline float getRightDistance() const;
  inline void setLeftDistance(float);
  inline void setRightDistance(float);
  virtual void postorderTraversal(TreeVisitor &);
  virtual void preorderTraversal(TreeVisitor &);
  virtual void recomputeNumLeaves();
  virtual bool reroot(const String &rootName);
  PhylogenyNode *getRemainingChild(float &getBranchLengthToo);
  inline void clearResolvedBit();
  inline void setResolvedBit();
  inline bool isResolved() const;
  virtual PhylogenyNode *clone();
  WhichChild whichChildIsThis(PhylogenyNode *);
  float getDistanceToChild(PhylogenyNode*);
  virtual PhylogenyNode *getChild(int i) {return i==0 ? left : right;}
};



/************************************************************************
                            class Phylogeny
*************************************************************************/
class Phylogeny
{
  PhylogenyNode *root;
  int numNodes;
  Array2D<float> pairwiseDistances; // between leaf taxa only

  void load(const String &filename);
  void pairwiseRecurs(PhylogenyNode *);
public:
  Phylogeny(PhylogenyNode *);
  Phylogeny(const String &filename);
  void save(const String &filename);
  PhylogenyNode *getRoot();
  bool reroot(const String &rootName);
  Phylogeny *clone();
  void printOn(ostream &);
  void printNewick(ostream &);
  void postorderTraversal(TreeVisitor &);
  void preorderTraversal(TreeVisitor &);
  void assignNodeIDs(); // also counts nodes
  void attachAlignment(MultSeqAlignment &);
  int getNumNodes() {return numNodes;}
  int getNumLeaves() {return root->getNumLeaves();}
  Vector<PhylogenyNode*> *gatherNodes(TraversalOrder=POSTORDER);
  void gatherNodes(Vector<PhylogenyNode*> &,TraversalOrder=POSTORDER);
  void gatherCladeMembers();
  Vector<PhylogenyBranch> *gatherBranches(); // DEPRECATED
  void gatherBranches(Vector<PhylogenyBranch> &); // DEPRECATED
  void collectBranches(Vector<PhylogenyBranch*> &,TraversalOrder=POSTORDER);
  void scaleBranches(double factor);
  float distanceBetween(PhylogenyNode *,PhylogenyNode *);
  float distanceBetweenLeaves(int ID1,int ID2);
  void computePairwiseDistances();
  float getSpannedDistance(int ID1,int ID2,int ID3);
  void prune(const Set<String> &keepTaxa);
  void constructBranches();
  void getPath(PhylogenyNode *from,PhylogenyNode *to,
	       Vector<PhylogenyBranch*> &);
  void getPath(PhylogenyNode *from,PhylogenyNode *to,Vector<PhylogenyNode*> &);
    // ^^^ includes beginning and end nodes
  PhylogenyNode *findNode(const String &name);
};

ostream &operator<<(ostream &,Phylogeny &);


/************************************************************************
                           class TreeVisitor
*************************************************************************/
class TreeVisitor {
public:
  virtual void processNode(InternalNode &) {}
  virtual void processNode(LeafNode &) {}
  virtual void processNode(RootNode &) {}
};



/************************************************************************
                           class NodeNumberer
*************************************************************************/
class NodeNumberer : public TreeVisitor
{
    int n;
public:
    NodeNumberer() : n(0) {}
    virtual void processNode(InternalNode &v) {v.setID(n++);}
    virtual void processNode(LeafNode &v) {v.setID(n++);}
    virtual void processNode(RootNode &v) {v.setID(n++);}
    int getNumNodes() {return n;}
};



/************************************************************************
                           class NodeCounter
*************************************************************************/
class NodeCounter : public TreeVisitor
{
    int n;
public:
    NodeCounter() : n(0) {}
    virtual void processNode(InternalNode &v) {n++;}
    virtual void processNode(LeafNode &v) {n++;}
    virtual void processNode(RootNode &v) {n++;}
    int getNumNodes() {return n;}
};



/************************************************************************
                        class AlignmentAttacher
*************************************************************************/
class AlignmentAttacher : public TreeVisitor
{
    MultSeqAlignment &A;
    int nextInternalNodeID, largestID;
public:
    AlignmentAttacher(MultSeqAlignment &);
    virtual void processNode(InternalNode &);
    virtual void processNode(LeafNode &);
    virtual void processNode(RootNode &);
};



/************************************************************************
                       class PhylogenyPruner
*************************************************************************/
class PhylogenyPruner : public TreeVisitor {
  const Set<String> &keepTaxa;
  void processChild(PhylogenyNode *&child,float &distance);
public:
  PhylogenyPruner(const Set<String> &keepTaxa);
  virtual void processNode(InternalNode &);
  virtual void processNode(LeafNode &);
  virtual void processNode(RootNode &);
};



/************************************************************************
                           inline methods
*************************************************************************/
ResidueCounter &PhylogenyNode::getCounter()
{
  return residueCounter;
}



PhylogenyNode *PhylogenyNode::getParent() const
{
  return parent;
}



PhylogenyNodeType PhylogenyNode::getNodeType() const
{
  return nodeType;
}



IndelHistory &PhylogenyNode::getIndelHistory()
{
    return indelHistory;
}



PhylogenyDecoration *&PhylogenyNode::getDecoration()
{
  return decoration;
}



void InternalNode::clearResolvedBit()
{
  resolvedBit=false;
}



void InternalNode::setResolvedBit()
{
  resolvedBit=true;
}



bool InternalNode::isResolved() const
{
  return resolvedBit;
}



void InternalNode::setLeftDistance(float d)
{
  leftDistance=d;
}



void InternalNode::setRightDistance(float d)
{
  rightDistance=d;
}



PhylogenyNode *InternalNode::getLeft()
{
  return left;
}



PhylogenyNode *InternalNode::getRight()
{
  return right;
}



void InternalNode::setLeft(PhylogenyNode *l)
{
  left=l;
}



void InternalNode::setRight(PhylogenyNode *r)
{
  right=r;
}



float InternalNode::getLeftDistance() const
{
  return leftDistance;
}



float InternalNode::getRightDistance() const
{
  return rightDistance;
}



WhichChild otherChild(WhichChild c)
{
  return c==LEFT ? RIGHT : LEFT;
}


#endif
