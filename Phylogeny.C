/****************************************************************
 Phylogeny.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "Phylogeny.H"
#include <iostream>
#include <fstream>
#include "BOOM/Exceptions.H"
using namespace std;
using namespace BOOM;


ostream &operator<<(ostream &os,PhylogenyNodeType t)
{
  switch(t) 
    {
    case ROOT_NODE:      os<<"ROOT_NODE"; break;
    case LEAF_NODE:      os<<"LEAF_NODE"; break;
    case INTERNAL_NODE:  os<<"INTERNAL_NODE"; break;
    }
  return os;
}


/************************************************************************
                            PhylogenyNode methods
*************************************************************************/
PhylogenyNode::PhylogenyNode(int numLeaves,PhylogenyNodeType nodeType) 
  : numLeaves(numLeaves),
    parent(NULL),
    nodeType(nodeType),
    decoration(NULL)
{
  // ctor
}



PhylogenyNode::PhylogenyNode(const PhylogenyNode &other)
  : numLeaves(other.numLeaves),
    nodeType(other.nodeType),
    parent(NULL),
    decoration(NULL)
{
  // copy ctor
}



void PhylogenyNode::setName(const String &n)
{
  name=n;
}



String &PhylogenyNode::getName()
{
  return name;
}



void PhylogenyNode::setID(int id)
{
  ID=id;
}



int PhylogenyNode::getID() const
{
  return ID;
}



void PhylogenyNode::setParent(PhylogenyNode *p)
{
  parent=p;
}



int PhylogenyNode::getNumLeaves() const
{
  return numLeaves;
}



void PhylogenyNode::setNumLeaves(int n)
{
  numLeaves=n;
}



PhylogenyNode *PhylogenyNode::load(istream &is)
{
  String line;
  line.getline(is);
  if(line.isEmpty() && is.eof()) throw "Error parsing phylogeny file";
  if(line=="LEAF") return new LeafNode(is);
  if(line=="INTERNAL") return new InternalNode(is);
  if(line=="ROOT") return new RootNode(is);
  throw "Error parsing phylogeny file";
}



WhichChild PhylogenyNode::whichChild()
{
  if(!parent) return LEFT;
  switch(parent->getNodeType()) 
    {
    case ROOT_NODE: return LEFT;
    case INTERNAL_NODE: {
      return static_cast<InternalNode*>(parent)->whichChildIsThis(this);
    }
    case LEAF_NODE: 
    default: 
      INTERNAL_ERROR;
    }
}



float PhylogenyNode::getDistanceToParent() const
{
  if(!parent) return 0.0;
  switch(parent->nodeType)
    {
    case ROOT_NODE: 
      return static_cast<RootNode*>(parent)->getBranchLength();
    case INTERNAL_NODE:
      return static_cast<InternalNode*>(parent)->getDistanceToChild(this);
    case LEAF_NODE: 
    default:
      INTERNAL_ERROR;
    }
}



void PhylogenyNode::allocateBranches(int n)
{
  branches.resize(n);
}



PhylogenyBranch &PhylogenyNode::getBranch(int i)
{
  return branches[i];
}



/************************************************************************
                            RootNode methods
*************************************************************************/
RootNode::RootNode(LeafNode &leaf,
		   float branchLength)
  : PhylogenyNode(0,ROOT_NODE),
    branchLength(branchLength),
    child(NULL),
    Pt(NULL)
{
  // ctor
}



RootNode::RootNode(float branchLength)
  : branchLength(branchLength),
    child(NULL),
    PhylogenyNode(0,ROOT_NODE),
    Pt(NULL)
{
  // ctor
}



RootNode::RootNode(istream &is)
  : PhylogenyNode(0,ROOT_NODE),
    Pt(NULL)
{
  // ctor

  String line;
  line.getline(is);
  setName(line);
  line.getline(is);
  branchLength=line.asFloat();
  child=PhylogenyNode::load(is);
  child->setParent(this);

  setNumLeaves(child->getNumLeaves());
}



RootNode::~RootNode()
{
  delete Pt;
}



float RootNode::getBranchLength() const
{
  return branchLength;
}



void RootNode::setBranchLength(float L)
{
  branchLength=L;
}



void RootNode::setSubstMatrix(NthOrdSubstMatrix *M)
{
  delete Pt;
  Pt=M;
}



NthOrdSubstMatrix *RootNode::getSubstMatrix()
{
  return Pt;
}



void RootNode::recomputeNumLeaves()
{
  child->recomputeNumLeaves();
  setNumLeaves(child->getNumLeaves());
}



void RootNode::output(ostream &os)
{
  os<<"ROOT"<<endl;
  os<<getName()<<endl;
  os<<branchLength<<endl;
  child->output(os);
}



void RootNode::printNewick(ostream &os)
{
  os<<"(";
  child->printNewick(os);
  os<<child->getName()<<":"<<branchLength<<")";
}



void RootNode::postorderTraversal(TreeVisitor &visitor)
{
  child->postorderTraversal(visitor);
  visitor.processNode(*this);
}



void RootNode::preorderTraversal(TreeVisitor &visitor)
{
  visitor.processNode(*this);
  child->preorderTraversal(visitor);
}



PhylogenyNode *RootNode::getChild()
{
  return child;
}



void RootNode::setChild(PhylogenyNode *child)
{
  this->child=child;
}



PhylogenyNode *RootNode::clone()
{
  RootNode *other=new RootNode(branchLength);
  other->child=child ? child->clone() : NULL;
  if(other->child) other->child->setParent(other);
  other->numLeaves=numLeaves;
  other->name=name;
  other->ID=ID;
  other->decoration=decoration ? decoration->clone() : NULL;
  other->cladeMembers=cladeMembers;
  other->indelHistory=indelHistory;
  return other;
}



/************************************************************************
                            LeafNode methods
*************************************************************************/
LeafNode::LeafNode() 
  : PhylogenyNode(1,LEAF_NODE) 
{
  // ctor
}



LeafNode::LeafNode(istream &is)
  : PhylogenyNode(1,LEAF_NODE)
{
  // ctor

  String name;
  name.getline(is);
  setName(name);
}



void LeafNode::recomputeNumLeaves()
{
  setNumLeaves(1);
}



void LeafNode::output(ostream &os)
{
  os<<"LEAF"<<endl;
  os<<getName()<<endl;
}



void LeafNode::printNewick(ostream &os)
{
}



void LeafNode::postorderTraversal(TreeVisitor &visitor)
{
  visitor.processNode(*this);
}



void LeafNode::preorderTraversal(TreeVisitor &visitor)
{
  visitor.processNode(*this);
}



bool LeafNode::reroot(const String &rootName)
{
  return getName()==rootName;
}



PhylogenyNode *LeafNode::clone()
{
  LeafNode *other=new LeafNode();
  other->name=name;
  other->numLeaves=numLeaves;
  other->ID=ID;
  other->indelHistory=indelHistory;
  other->decoration=decoration ? decoration->clone() : NULL;
  other->cladeMembers=cladeMembers;
  return other;
}



/************************************************************************
                          InternalNode methods
*************************************************************************/
InternalNode::InternalNode(PhylogenyNode *left,PhylogenyNode *right) 
  : left(left), 
    right(right), 
    PhylogenyNode(left->getNumLeaves()+right->getNumLeaves(),INTERNAL_NODE),
    leftDistance(0),
    rightDistance(0),
    resolvedBit(false),
    leftMatrix(NULL),
    rightMatrix(NULL)
{
  // ctor

  left->setParent(this);
  right->setParent(this);
}



InternalNode::InternalNode(istream &is)
  : PhylogenyNode(0,INTERNAL_NODE),
    leftMatrix(NULL),
    rightMatrix(NULL)
{
  // ctor

  String line;
  line.getline(is);
  setName(line);
  line.getline(is);
  leftDistance=line.asFloat();
  left=PhylogenyNode::load(is);
  left->setParent(this);
  line.getline(is);
  rightDistance=line.asFloat();
  right=PhylogenyNode::load(is);
  right->setParent(this);

  setNumLeaves(left->getNumLeaves()+right->getNumLeaves());
}



InternalNode::~InternalNode()
{
  delete leftMatrix;
  delete rightMatrix;
}



void InternalNode::setLeftSubstMatrix(NthOrdSubstMatrix *M)
{
  delete leftMatrix;
  leftMatrix=M;
}



void InternalNode::setRightSubstMatrix(NthOrdSubstMatrix *M)
{
  delete rightMatrix;
  rightMatrix=M;
}



NthOrdSubstMatrix *InternalNode::getLeftSubstMatrix()
{
  return leftMatrix;
}



NthOrdSubstMatrix *InternalNode::getRightSubstMatrix()
{
  return rightMatrix;
}



void InternalNode::recomputeNumLeaves()
{
  left->recomputeNumLeaves();
  right->recomputeNumLeaves();
  setNumLeaves(left->getNumLeaves()+right->getNumLeaves());
}



void InternalNode::output(ostream &os)
{
  os<<"INTERNAL"<<endl;
  os<<getName()<<endl;
  os<<leftDistance<<endl;
  left->output(os);
  os<<rightDistance<<endl;
  right->output(os);
}



void InternalNode::printNewick(ostream &os)
{
  os<<"(";
  left->printNewick(os);
  os<<left->getName()<<":"<<leftDistance<<",";
  right->printNewick(os);
  os<<right->getName()<<":"<<rightDistance<<")";
}



void InternalNode::postorderTraversal(TreeVisitor &visitor)
{
  left->postorderTraversal(visitor);
  right->postorderTraversal(visitor);
  visitor.processNode(*this);
}



void InternalNode::preorderTraversal(TreeVisitor &visitor)
{
  visitor.processNode(*this);
  left->preorderTraversal(visitor);
  right->preorderTraversal(visitor);
}



bool InternalNode::reroot(const String &rootName)
{
  if(left->reroot(rootName))
    {
      reroot(left,LEFT,leftDistance);
      return true;
    }
  if(right->reroot(rootName))
    {
      reroot(right,RIGHT,rightDistance);
      return true;
    }
  return false;
}



WhichChild InternalNode::whichChildIsThis(PhylogenyNode *child)
{
  return left==child ? LEFT : RIGHT;
}



PhylogenyNode *InternalNode::getOtherChild(WhichChild c)
{
  return c==LEFT ? right : left;
}



float &InternalNode::getDistance(WhichChild c)
{
  return c==LEFT ? leftDistance : rightDistance;
}



void InternalNode::reroot(PhylogenyNode *child,WhichChild whichChild,
			  float branchLength)
{
  // Get the node that was my parent, and the branch length to it
  PhylogenyNode *parent=getParent();
  bool wasRoot=!parent;
  float distanceToMe=0;
  if(parent) 
    {
      InternalNode *node=(InternalNode*) parent;
      WhichChild whichChildAmI=node->whichChildIsThis(this);
      distanceToMe=node->getDistance(whichChildAmI);
    }

  // Plug my old parent in as my new child
  if(!wasRoot)
    if(whichChild==LEFT) 
      { left=parent; leftDistance=distanceToMe; }
    else 
      { right=parent; rightDistance=distanceToMe; }
  
  // Make my old child become my new parent
  if(child->getNodeType()==LEAF_NODE)
    {
      LeafNode *leaf=(LeafNode*) child;
      RootNode *root=new RootNode(*leaf,branchLength);
      root->setName(child->getName());
      delete child;
      root->setChild(this);
      setParent(root);
      child=root;
    }
  else setParent(child);

  if(wasRoot)
    {
      // This was previously the root of the tree; it will now have only
      // one child, which is illegal in a binary tree, so we need to
      // re-route all links around it so it can be deleted.

      // First, tell my new parent that I am no longer its left/right child:
      PhylogenyNode *myOtherChild=getOtherChild(whichChild);
      WhichChild otherChild=(whichChild==LEFT ? RIGHT : LEFT);
      double otherDistance=(otherChild==LEFT ? leftDistance : rightDistance);
      switch(child->getNodeType())
	{
	case INTERNAL_NODE:
	  {
	    InternalNode *node=(InternalNode*) child;
	    if(node->whichChildIsThis(this)==LEFT)
	      {
		node->left=myOtherChild;
		node->leftDistance+=otherDistance;//###leftDistance;
	      }
	    else
	      {
		node->right=myOtherChild;
		node->rightDistance+=otherDistance;//###rightDistance;
	      }
	  }
	  break;
	case ROOT_NODE:
	  {
	    RootNode *node=(RootNode*)child;
	    node->setChild(myOtherChild);
	    node->setBranchLength(node->getBranchLength()+otherDistance);
	  }
	  break;
	}

      // Now tell my remaining child who his new parent is
      myOtherChild->setParent(child);

      // Now this node can be safely deleted in Phylogeny::reroot().
    }
}



PhylogenyNode *InternalNode::getRemainingChild(float &getBranchLengthToo)
{
  if(left)
    {
      getBranchLengthToo=leftDistance;
      return left;
    }
  getBranchLengthToo=rightDistance;
  return right;
}



PhylogenyNode *InternalNode::clone()
{
  InternalNode *other=new InternalNode(left->clone(),right->clone());
  other->leftDistance=leftDistance;
  other->rightDistance=rightDistance;
  if(other->left) other->left->setParent(other);
  if(other->right) other->right->setParent(other);
  other->numLeaves=numLeaves;
  other->name=name;
  other->ID=ID;
  other->indelHistory=indelHistory;
  other->decoration=decoration ? decoration->clone() : NULL;
  other->cladeMembers=cladeMembers;
  return other;
}



float InternalNode::getDistanceToChild(PhylogenyNode *child)
{
  return getDistance(whichChildIsThis(child));
}



int InternalNode::numChildren()
{
  if(left) 
    if(right) return 2;
    else return 1;
  else
    if(right) return 1;
    else return 0;
}



/************************************************************************
                           Phylogeny methods
*************************************************************************/
Phylogeny::Phylogeny(PhylogenyNode *root)
  : root(root), numNodes(0)
{
  // ctor
}



Phylogeny::Phylogeny(const String &filename)
  : numNodes(0)
{
  // ctor

  load(filename);
}



PhylogenyNode *Phylogeny::getRoot()
{
  return root;
}



void Phylogeny::load(const String &filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw filename+" could not be loaded in Phylogeny::load()";
  root=PhylogenyNode::load(is);
  assignNodeIDs();
}



bool Phylogeny::reroot(const String &rootName)
{
  PhylogenyNode *oldRoot=root;
  if(!root->reroot(rootName)) return false;
  while(1)
    {
      PhylogenyNode *parent=root->getParent();
      if(!parent) break;
      root=parent;
    }
  if(root!=oldRoot && oldRoot->getNodeType()==INTERNAL_NODE) {
    delete oldRoot;
  }
  root->recomputeNumLeaves();
  return true;
}



Phylogeny *Phylogeny::clone()
{
  Phylogeny *other=new Phylogeny(root->clone());
  other->root->recomputeNumLeaves();
  NodeCounter counter;
  other->postorderTraversal(counter);
  other->numNodes=counter.getNumNodes();
  return other;
}



void Phylogeny::postorderTraversal(TreeVisitor &visitor)
{
  root->postorderTraversal(visitor);
}



void Phylogeny::preorderTraversal(TreeVisitor &visitor)
{
    root->preorderTraversal(visitor);
}



void Phylogeny::printOn(ostream &os)
{
  root->output(os);
}



void Phylogeny::printNewick(ostream &os)
{
  root->printNewick(os);
  os<<root->getName();
}



ostream &operator<<(ostream &os,Phylogeny &phy)
{
  phy.printOn(os);
  return os;
}



void Phylogeny::assignNodeIDs()
{
  NodeNumberer numberer;
  postorderTraversal(numberer);
  numNodes=numberer.getNumNodes();
}



void Phylogeny::attachAlignment(MultSeqAlignment &alignment)
{
  AlignmentAttacher A(alignment);
  postorderTraversal(A);
}



Vector<PhylogenyNode*> *Phylogeny::gatherNodes(TraversalOrder order)
{
  Vector<PhylogenyNode*> *nodes=new Vector<PhylogenyNode*>;
  gatherNodes(*nodes,order);
  return nodes;
}



void Phylogeny::gatherNodes(Vector<PhylogenyNode*> &nodes,
			    TraversalOrder order)
{
  nodes.clear();
  struct Gatherer : public TreeVisitor {
    void processNode(InternalNode &node) {nodes.push_back(&node);}
    void processNode(LeafNode &node) {nodes.push_back(&node);}
    void processNode(RootNode &node) {nodes.push_back(&node);}
    Gatherer(Vector<PhylogenyNode*> &nodes) : nodes(nodes) {}
    Vector<PhylogenyNode*> &nodes;
  } gatherer(nodes);
  switch(order) 
    {
    case POSTORDER: postorderTraversal(gatherer); break;
    case PREORDER: preorderTraversal(gatherer); break;
    case INORDER: throw "option INORDER not currently supported in Phylogeny";
    }
}



Vector<PhylogenyBranch> *Phylogeny::gatherBranches()
{
  Vector<PhylogenyBranch> *v=new Vector<PhylogenyBranch>;
  gatherBranches(*v);
  return v;
}



void Phylogeny::collectBranches(Vector<PhylogenyBranch*> &branches,
				TraversalOrder order)
{
  struct Gatherer : public TreeVisitor {
    Vector<PhylogenyBranch*> &branches;
    Gatherer(Vector<PhylogenyBranch*> &branches) : branches(branches) {}
    void processNode(InternalNode &node) 
    {
      branches.push_back(&node.getBranch(0));
      branches.push_back(&node.getBranch(1));
    }
    void processNode(RootNode &node) 
    {
      branches.push_back(&node.getBranch(0));
    }
  } gatherer(branches);
  switch(order) 
    {
    case PREORDER:   preorderTraversal(gatherer);  break;
    case POSTORDER:  postorderTraversal(gatherer);  break;
    case INORDER:    INTERNAL_ERROR;
    }
}



void Phylogeny::gatherBranches(Vector<PhylogenyBranch> &branches)
{
  struct Gatherer : public TreeVisitor {
    void processNode(InternalNode &node) 
    {
      branches.push_back(PhylogenyBranch(&node,node.getLeft(),
					 branches.size(),
					 node.getLeftDistance()));
      branches.push_back(PhylogenyBranch(&node,node.getRight(),
					 branches.size(),
					 node.getRightDistance())); 
    }
    void processNode(RootNode &node) 
    { 
      branches.push_back(PhylogenyBranch(&node,node.getChild(),
					 branches.size(),
					 node.getBranchLength()));
    }
    void processNode(LeafNode &node) {}
    Gatherer(Vector<PhylogenyBranch> &branches) : branches(branches) {}
    Vector<PhylogenyBranch> &branches;
  } gatherer(branches);
  postorderTraversal(gatherer);
}




//void gatherBranches(Vector<PhylogenyBranch> &);
void Phylogeny::scaleBranches(double factor)
{
  Vector<PhylogenyBranch> branches;
  gatherBranches(branches);
  Vector<PhylogenyBranch>::iterator cur=branches.begin(), end=branches.end();
  for(; cur!=end ; ++cur) {
    PhylogenyBranch &branch=*cur;
    double oldLen=branch.getLength();
    double newLen=oldLen*factor;
    branch.changeLength(newLen);
  }
}



float Phylogeny::distanceBetween(PhylogenyNode *a,PhylogenyNode *b)
{
  // First, construct path to root, with annotated distance sums
  Map<PhylogenyNode*,float> distances;
  float distance=0.0;
  PhylogenyNode *A=a;
  while(A) {
    distances[A]=distance;
    distance+=A->getDistanceToParent();
    A=A->getParent();
    if(A==b) return distance;
  }
  
  // Now search upward rom the other node, till we come to the path
  // connecting the first node to the root
  distance=0.0;
  while(!distances.isDefined(b)) {
    distance+=b->getDistanceToParent();
    b=b->getParent();
  }
  return distance+distances[b];
}




void Phylogeny::prune(const Set<String> &keepTaxa)
{
  /*
  cout<<"keeping: "<<endl;
  Set<String>::iterator cur=keepTaxa.begin(), end=keepTaxa.end();
  for(; cur!=end ; ++cur) cout<<*cur<<endl;
  cout<<"."<<endl;
  */

  PhylogenyPruner pruner(keepTaxa);
  root->postorderTraversal(pruner);
  if(root->getNodeType()==INTERNAL_NODE) {
    InternalNode *node=static_cast<InternalNode*>(root);
    switch(node->numChildren()) {
    case 0: INTERNAL_ERROR;
    case 1: {
      float dummy;
      root=node->getRemainingChild(dummy);
      root->setParent(NULL);
      delete node;
    }
    case 2: break; // do nothing
    }
  }
  NodeCounter counter;
  root->postorderTraversal(counter);
  numNodes=counter.getNumNodes();
}



void Phylogeny::save(const String &filename)
{
  ofstream os(filename.c_str());
  if(!os.good()) throw String("Error writing to file: ")+filename;
  os<<*this<<endl;
}



void Phylogeny::computePairwiseDistances()
{
  gatherCladeMembers();
  pairwiseDistances.resize(numNodes,numNodes);
  pairwiseDistances.setAllTo(0.0);

  struct Visitor {
    BitSet accum;
    Array2D<float> &dist;
    int numTaxa;
    Visitor(Array2D<float> &dist) 
      : accum(dist.getFirstDim()), dist(dist) {
      numTaxa=dist.getFirstDim();
    }
    void add(float d,BitSet who) {
      for(int i=0 ; i<numTaxa ; ++i)
	if(who.isMember(i)) dist[i][i]+=d;
    }
    void recurs(PhylogenyNode *node) {
      switch(node->getNodeType())
	{
	case ROOT_NODE: {
	  RootNode *root=static_cast<RootNode*>(node);
	  accum.addMember(root->getID());
	  recurs(root->getChild());
	  }
	  break;
	case LEAF_NODE: {
	  LeafNode *leaf=static_cast<LeafNode*>(node);
	  int ID=leaf->getID();
	  for(int i=0 ; i<numTaxa ; ++i)
	    if(accum.isMember(i)) 
	      dist[ID][i]=dist[i][ID]=dist[i][i];
	  accum.addMember(ID);
	  }
	  break;
	case INTERNAL_NODE: {
	  InternalNode *intNode=static_cast<InternalNode*>(node);
	  PhylogenyNode *left=intNode->getLeft(), *right=intNode->getRight();
	  float ld=intNode->getLeftDistance(), rd=intNode->getRightDistance();
	  accum.addMember(intNode->getID());
	  BitSet parental=accum;
	  add(ld,parental);
	  recurs(left);
	  add(-ld,parental);
	  add(ld,left->getCladeMembers());
	  add(rd,accum);
	  BitSet nonRight=accum;
	  recurs(right);
	  add(-rd,nonRight);
	  add(rd,right->getCladeMembers());
	  }
	  break;
	}
    }
  } visitor(pairwiseDistances);
  visitor.recurs(root);
  for(int i=0 ; i<numNodes ; ++i) pairwiseDistances[i][i]=0;
}



void Phylogeny::gatherCladeMembers()
{
  struct TV : public TreeVisitor {
    int numTaxa;
    TV(int numTaxa) : numTaxa(numTaxa) {}
    void processNode(InternalNode &node) {
      int id=node.getID();
      BitSet &M=node.getCladeMembers();
      M=node.getLeft()->getCladeMembers();
      M+=node.getRight()->getCladeMembers();
      M.addMember(id);
    }
    void processNode(LeafNode &node) {
      int id=node.getID();
      BitSet &M=node.getCladeMembers();
      M.setSize(numTaxa);
      M.addMember(id);
    }
    void processNode(RootNode &node) {
      int id=node.getID();
      BitSet &M=node.getCladeMembers();
      M=node.getChild()->getCladeMembers();
      M.addMember(id);
    }
  } visitor(numNodes);
  postorderTraversal(visitor);
}



float Phylogeny::getSpannedDistance(int ID1,int ID2,int ID3)
{
  if(pairwiseDistances.getFirstDim()==0) computePairwiseDistances();
  return 
    (pairwiseDistances[ID1][ID2]+
     pairwiseDistances[ID2][ID3]+
     pairwiseDistances[ID3][ID1])/2;
}



float Phylogeny::distanceBetweenLeaves(int ID1,int ID2)
{
  if(pairwiseDistances.getFirstDim()==0) computePairwiseDistances();
  return pairwiseDistances[ID1][ID2];
}



void Phylogeny::constructBranches()
{
  struct TV : public TreeVisitor {
    int nextID;
    TV() : nextID(0) {}
    void processNode(RootNode &node) {
      node.allocateBranches(1);
      PhylogenyBranch *branch=&node.getBranch(0);
      branch->setParent(&node);
      branch->setChild(node.getChild());
      branch->changeLength(node.getBranchLength());
      branch->setBranchID(nextID++);
    }
    void processNode(InternalNode &node) {
      node.allocateBranches(2);
      PhylogenyBranch *leftBranch=&node.getBranch(0);
      PhylogenyBranch *rightBranch=&node.getBranch(1);
      leftBranch->setParent(&node);
      rightBranch->setParent(&node);
      leftBranch->setChild(node.getLeft());
      rightBranch->setChild(node.getRight());
      leftBranch->changeLength(node.getLeftDistance());
      rightBranch->changeLength(node.getRightDistance());
      leftBranch->setBranchID(nextID++);
      rightBranch->setBranchID(nextID++);
    }
  } visitor;
  preorderTraversal(visitor);
}



void PhylogenyNode::getPathToRoot(Vector<PhylogenyNode*> &path)
{
  path.push_back(this);
  PhylogenyNode *node=this;
  while(node->parent) {
    node=node->parent;
    path.push_back(node);
  }
}



void Phylogeny::getPath(PhylogenyNode *from,PhylogenyNode *to,
			Vector<PhylogenyNode*> &path)
{
  Vector<PhylogenyNode*> fromPath, toPath;
  from->getPathToRoot(fromPath);
  to->getPathToRoot(toPath);
  PhylogenyNode *commonAncestor=fromPath.pop();
  toPath.pop();
  while(fromPath.back()==toPath.back()) {
    commonAncestor=fromPath.pop();
    toPath.pop();
  }
  int n1=fromPath.size(), n2=toPath.size();
  for(int i=0 ; i<n1 ; ++i) path.push_back(fromPath[i]);
  path.push_back(commonAncestor);
  for(int i=n2-1 ; i>=0 ; --i) path.push_back(toPath[i]);
}



void Phylogeny::getPath(PhylogenyNode *from,PhylogenyNode *to,
			Vector<PhylogenyBranch*> &path)
{
  path.clear();
  Vector<PhylogenyNode*> nodes;
  getPath(from,to,nodes);
  int n=nodes.size();
  for(int i=0 ; i<n-1 ; ++i) {
    PhylogenyNode *thisNode=nodes[i], *nextNode=nodes[i+1];
    if(thisNode->getParent()==nextNode) {
      PhylogenyNode *temp=thisNode;
      thisNode=nextNode;
      nextNode=temp;
    }
    PhylogenyBranch *branch=&thisNode->getBranch(nextNode->whichChild());
    if(branch->getChild()!=nextNode) INTERNAL_ERROR;
    path.push_back(branch);
  }
}



PhylogenyNode *Phylogeny::findNode(const String &name)
{
  struct Finder : public TreeVisitor {
    PhylogenyNode *found;
    String name;
    Finder(const String &name) : name(name), found(NULL) {}
    void processNode(InternalNode &v) {process(v);}
    void processNode(LeafNode &v) {process(v);}
    void processNode(RootNode &v) {process(v);}
    void process(PhylogenyNode &v) {if(v.getName()==name) found=&v;}
  } finder(name);
  postorderTraversal(finder);
  return finder.found;
}



/************************************************************************
                        AlignmentAttacher methods
*************************************************************************/
AlignmentAttacher::AlignmentAttacher(MultSeqAlignment &A) 
    : A(A), nextInternalNodeID(A.getNumTracks()), largestID(0) 
{
}



void AlignmentAttacher::processNode(InternalNode &V) 
{
    AlignmentSeq &track=A.findOrCreateTrack(V.getName());
    V.setID(track.getID());
}



void AlignmentAttacher::processNode(LeafNode &V)
{
  //int ID=A.getTrackByName(V.getName()).getID();
  int ID=A.findOrCreateTrack(V.getName()).getID();
  if(ID>largestID) largestID=ID;
  V.setID(ID);
  A.getIthTrack(ID).extendToLength(A.getLength(),A.getGapSymbol());
}



void AlignmentAttacher::processNode(RootNode &V)
{
  //int ID=A.getTrackByName(V.getName()).getID();
    int ID=A.findOrCreateTrack(V.getName()).getID();
    if(ID>largestID) largestID=ID;
    V.setID(ID);
}



/************************************************************************
                         PhylogenyBranch methods
*************************************************************************/

PhylogenyBranch::PhylogenyBranch(PhylogenyNode *parent,PhylogenyNode *child,
				 int branchID,double length)
  : parent(parent), child(child), branchID(branchID), length(length),
    decoration(NULL)
{
  // ctor
}



PhylogenyBranch::PhylogenyBranch()
  : parent(NULL), child(NULL), branchID(-1), length(-1.0), decoration(NULL)
{
  // default ctor
}



PhylogenyNode *PhylogenyBranch::getParent() const
{
  return parent;
}



PhylogenyNode *PhylogenyBranch::getChild() const
{
  return child;
}



int PhylogenyBranch::getBranchID() const
{
  return branchID;
}



double PhylogenyBranch::getLength() const
{
  return length;
}



void PhylogenyBranch::changeLength(double newLength)
{
  length=newLength;
  switch(parent->getNodeType())
    {
    case ROOT_NODE:
      static_cast<RootNode*>(parent)->setBranchLength(newLength);
      break;
    case LEAF_NODE:
      INTERNAL_ERROR;
      break;
    case INTERNAL_NODE:
      {
	InternalNode *iParent=static_cast<InternalNode*>(parent);
	WhichChild whichChild=iParent->whichChildIsThis(child);
	if(whichChild==LEFT) iParent->setLeftDistance(newLength);
	else iParent->setRightDistance(newLength);
      }
      break;
    }
}


BranchDecoration *&PhylogenyBranch::getDecoration()
{
  return decoration;
}



/************************************************************************
                       PhylogenyPruner methods
*************************************************************************/

PhylogenyPruner::PhylogenyPruner(const Set<String> &keepTaxa)
  : keepTaxa(keepTaxa)
{
  // ctor
}



void PhylogenyPruner::processNode(InternalNode &node)
{
  PhylogenyNode *left=node.getLeft(), *right=node.getRight();
  float leftDistance=node.getLeftDistance();
  float rightDistance=node.getRightDistance();
  processChild(left,leftDistance);
  processChild(right,rightDistance);
  node.setLeft(left);
  node.setRight(right);
  node.setLeftDistance(leftDistance);
  node.setRightDistance(rightDistance);
  if(left) left->setParent(&node);
  if(right) right->setParent(&node);
}



void PhylogenyPruner::processChild(PhylogenyNode *&child,float &distance)
{
  if(child->getNodeType()==LEAF_NODE) {
    if(!keepTaxa.isMember(child->getName())) {
      delete child;
      child=NULL;
    }
  }
  else {
    InternalNode *node=static_cast<InternalNode*>(child);
    switch(node->numChildren()) {
    case 0:
      delete child;
      child=NULL;
      break;
    case 1:
      {
	float additionalDistance;
	PhylogenyNode *grandchild=
	  node->getRemainingChild(additionalDistance);
	PhylogenyNode *parent=node->getParent();
	delete node;
	child=grandchild;
	distance+=additionalDistance;
      }
      break;
    case 2: break; // do nothing
    }
  }
}



void PhylogenyPruner::processNode(LeafNode &node)
{
  // nothing to do
}



void PhylogenyPruner::processNode(RootNode &root)
{
  PhylogenyNode *child=root.getChild();
  float branchLen=root.getBranchLength();
  processChild(child,branchLen);
  root.setChild(child);
  root.setBranchLength(branchLen);
  if(child) child->setParent(&root);
}



