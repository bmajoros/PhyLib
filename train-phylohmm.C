/****************************************************************
 train-phylohmm.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/


THIS IS OBSOLETE


#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Vector.H"
#include "BOOM/Map.H"
#include "BOOM/Random.H"
#include "BOOM/Constants.H"
#include "BOOM/ConfigFile.H"
#include "BOOM/GSL/Optimizer.H"
#include "BOOM/Environment.H"
#include "BOOM/MolecularSequenceType.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/AminoAlphabet.H"
#include "Phylogeny.H"
#include "NthOrdRateMatrix.H"
#include "RateMatrixType.H"
#include "NthOrdFelsenstein.H"
using namespace std;
using namespace BOOM;


enum GradientMethod {
  DIFFERENCING_METHOD,
  ANALYTICAL_METHOD
};

/****************************************************************
                     class ObjectiveFunction
 ****************************************************************/
class ObjectiveFunction : public GSL::ObjectiveFunction
{
  GradientMethod gradientMethod;
  MultSeqAlignment &alignment;
  Phylogeny &phylogeny;
  NthOrdRateMatrix &Q;
  int dimensionality;
  double worstLikelihood, likelihood;
  GSL::Optimizer **optimizer;
  double dx;
  bool treeIsFixed;
  int numThreads;
  MolecularSequenceType seqType;
  ContextType contextType;
  Map<const GSL::Vector,double> cache;
  int cacheHits;
  int numEvaluations;

  void differencingMethod(const GSL::Vector &currentPoint,
			    GSL::Vector &gradient);
  void analyticalMethod(const GSL::Vector &currentPoint,
			  GSL::Vector &gradient);
  bool pointIsValid(const GSL::Vector &point);
public:
  ObjectiveFunction(MultSeqAlignment &,Phylogeny &,NthOrdRateMatrix &,
		    GradientMethod,int dimensionality,
		    GSL::Optimizer **,double gradientEpsilon=1e-6,
		    bool treeIsFixed=false,int numThreads=1,
		    MolecularSequenceType=DNA,ContextType=CT_RCO);
  virtual double f(const GSL::Vector &currentPoint);
  virtual void gradient(const GSL::Vector &currentPoint,
			GSL::Vector &gradient);
  void copyIn(GSL::Vector &);
  void copyOut(const GSL::Vector &);
  int getNumEvaluations() const {return numEvaluations;}
  int getNumCacheHits() const {return cacheHits;}
};



/****************************************************************
                       class Application
 ****************************************************************/
class Application
{
  MultSeqAlignment *alignment;
  Phylogeny *phylogeny;
  RateMatrixType matrixType;
  NthOrdRateMatrix *Q;
  int dimensionality; // total number of parameters to optimize
  GSL::Optimizer *optimizer;
  String optimizerTypeStr, stoppingCriterionStr;
  double likelihoodThreshold, gradientEpsilon, stepSize;
  bool treeIsFixed;
  int numThreads;
  bool dual;
  Alphabet *alphabet;
  Symbol gapSymbol;
  MolecularSequenceType seqType;
  int order;
  float noiseSigma; // for generating random (higher-order) rate matrix
  ContextType contextType;

  void getEquilibriumFreqs(MultiAlignment &,
			   Array1D<double> &equilibriumFreqs);
  void optimize(GradientMethod,int maxIterations);
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
                    class GetBranchLengths
 ****************************************************************/
class GetBranchLengths : public TreeVisitor
{
  int nextIndex;
  GSL::Vector &v;
public:
  GetBranchLengths(GSL::Vector &);
  virtual void processNode(InternalNode &);
  virtual void processNode(LeafNode &) {}
  virtual void processNode(RootNode &);
  int getNextIndex() {return nextIndex;}
};



/****************************************************************
                    class SetBranchLengths
 ****************************************************************/
class SetBranchLengths : public TreeVisitor
{
  int nextIndex;
  const GSL::Vector &v;
public:
  SetBranchLengths(const GSL::Vector &);
  virtual void processNode(InternalNode &);
  virtual void processNode(LeafNode &) {}
  virtual void processNode(RootNode &);
  int getNextIndex() {return nextIndex;}
};



/****************************************************************
                      Application methods
 ****************************************************************/



Application::Application()
    : Q(NULL), alignment(NULL), phylogeny(NULL), optimizer(NULL),
      treeIsFixed(false), noiseSigma(0.05)
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"m:td");
    if(cmd.numArgs()!=5)
      throw String(
"\ntrain-phylohmm <in.phy> <in.maf> <out.phy> <out.matrix> <config-file>\n\
\n\
options: \n\
         -m matrix = seed the matrix rather than initializing randomly\n\
         -t = don't modify the tree parameters\n\
         -d = use dual contexts\n\
\n"
);
    String phyFile=cmd.arg(0);
    String mafFile=cmd.arg(1);
    String outPhy=cmd.arg(2);
    String outMatrix=cmd.arg(3);
    String configFilename=cmd.arg(4);
    treeIsFixed=cmd.option('t');
    randomize();
    numThreads=Environment::lookup("NUM_CPUS");
    if(numThreads<1) numThreads=1;
    dual=cmd.option('d');
    
    // Process the configuration file
    ConfigFile configFile(configFilename);
    String Mtype=configFile.lookupOrDie("matrix-type");
    int maxIterations=configFile.getIntOrDie("max-iterations");
    optimizerTypeStr=configFile.lookupOrDie("optimizer");
    stoppingCriterionStr=configFile.lookupOrDie("stopping-criterion");
    likelihoodThreshold=configFile.getDoubleOrDie("likelihood-threshold");
    gradientEpsilon=configFile.getDoubleOrDie("gradient-epsilon");
    stepSize=configFile.getDoubleOrDie("step-size");
    seqType=seqTypeFromString(configFile.lookupOrDie("sequence-type"));
    order=configFile.getIntOrDie("order");
    if(configFile.isDefined("noise-sigma")) 
      noiseSigma=configFile.getFloatOrDie("noise-sigma");
    contextType=contextTypeFromString(configFile.lookupOrDie("context-type"));

    // Load the alignments
    Vector<MultiAlignment*> alignments;
    MultiAlignment::loadMAF(mafFile,alignments);
    MultiAlignment *alignment=MultiAlignment::combine(alignments,true);
    alignment->toupper();
    alphabet= seqType==DNA ? 
      &(Alphabet&) DnaAlphabet::global : &(Alphabet&) AminoAlphabet::global;
    gapSymbol=seqType==DNA ? alphabet->lookup('N') : alphabet->lookup('*');
    MultSeqAlignment &A=*new MultSeqAlignment(*alignment,*alphabet,gapSymbol);
    this->alignment=&A;

    // Load the phylogeny
    phylogeny=new Phylogeny(phyFile);
    phylogeny->attachAlignment(A);

    // Get equilibrium frequencies
    Array1D<double> equilibriumFreqs(5);
    getEquilibriumFreqs(*alignment,equilibriumFreqs);

    // Load or generate the rate matrix
    if(cmd.option('m')) 
      {
	Q=NthOrdRateMatrix::load(cmd.optParm('m'));
	matrixType=Q->getMatrixType();
      }
    else 
      {
	int ord=dual ? 2*order : order;
	matrixType=rateMatrixTypeFromString(Mtype);
	Q=NthOrdRateMatrix::random(dual,noiseSigma,seqType,ord,
				   matrixType,equilibriumFreqs);
      }
    dimensionality=Q->numParameters()*Q->getNumMatrices();
    if(!treeIsFixed) dimensionality+=phylogeny->getNumNodes()-1;

    // Report initial likelihood, for comparison later
    //NthOrdFelsenstein fel(contextType,*phylogeny,Q,A,seqType,numThreads);
    //cout<<"Initial likelihood: "<<fel.logLikelihood()<<endl;

    // Optimize
    GradientMethod gradientMethod=DIFFERENCING_METHOD;
    optimize(gradientMethod,maxIterations);

    // Generate output
    ofstream os(outPhy.c_str());
    phylogeny->printOn(os);
    Q->save(outMatrix);

    cout<<"Done."<<endl;
    return 0;
  }



void Application::getEquilibriumFreqs(MultiAlignment &alignment,
				      Array1D<double> &equilibriumFreqs)
{
  int n=alignment.getNumTracks();
  int sampleSize=0;
  equilibriumFreqs.resize(5);
  equilibriumFreqs.setAllTo(0.0);
  for(int i=0 ; i<n ; ++i)
    {
      AlignmentTrack &track=alignment.getIthTrack(i);
      const String &seq=track.getSeq();
      int L=seq.length();
      for(int i=0 ; i<L ; ++i)
	switch(seq[i])
	  {
	  case 'A': ++equilibriumFreqs[0]; ++sampleSize; break;
	  case 'C': ++equilibriumFreqs[1]; ++sampleSize; break;
	  case 'G': ++equilibriumFreqs[2]; ++sampleSize; break;
	  case 'T': ++equilibriumFreqs[4]; ++sampleSize; break;
	  case 'N': ++equilibriumFreqs[3]; ++sampleSize; break;
	  }
    }
  for(int i=0 ; i<5 ; ++i)
    equilibriumFreqs[i]/=sampleSize;
}



void Application::optimize(GradientMethod gradientMethod,int maxIterations)
{
  ObjectiveFunction f(*alignment,*phylogeny,*Q,gradientMethod,
		      dimensionality,&optimizer,gradientEpsilon,
                      treeIsFixed,numThreads,seqType,contextType);
  GSL::Vector initialPoint, lowerBounds, upperBounds;
  initialPoint.resize(dimensionality);
  f.copyIn(initialPoint);
  cout<<"INITIAL POINT: "<<initialPoint<<endl;
  GSL::OptimizerType optimizerType=
      GSL::stringToOptimizerType(optimizerTypeStr);
  GSL::StoppingCriterion stoppingCriterion=
      GSL::stringToStoppingCriterion(stoppingCriterionStr);

  lowerBounds.resize(dimensionality);
  upperBounds.resize(dimensionality);
  lowerBounds.setAllTo(0.0001);
  upperBounds.setAllTo(0.99);
  int numMatrixParms=Q->numParameters()*Q->getNumMatrices();
  int numBranches=dimensionality-numMatrixParms;
  for(int i=0 ; i<numBranches ; ++i) upperBounds[i]=1000000;

  /*optimizer=new GSL::Optimizer(optimizerType,f,initialPoint,stepSize,
                               stoppingCriterion,likelihoodThreshold,
                               maxIterations);*/
  optimizer=new GSL::ConstrainedOptimizer(f,initialPoint,lowerBounds,
					  upperBounds,stepSize,
					  stoppingCriterion,
					  likelihoodThreshold,maxIterations);
  optimizer->run();
  cout<<"optimization complete."<<endl;
  const GSL::Vector &optimalPoint=optimizer->getOptimalPoint();
  //###delete optimizer;
  f.copyOut(optimalPoint);
  cout<<"OPTIMAL POINT: "<<optimalPoint<<endl;
  cout<<"logL increase: -"
      <<((GSL::ConstrainedOptimizer*)optimizer)->getInitialValue()
      <<" to -"
      <<((GSL::ConstrainedOptimizer*)optimizer)->getOptimalValue()
    //<<" in "<<optimizer->iterationsUsed()<<" iterations"
      <<" in "<<f.getNumEvaluations()<<" evals + "
      <<f.getNumCacheHits()<<" cache hits"
      <<endl;
}



/****************************************************************
                    ObjectiveFunction methods
 ****************************************************************/

ObjectiveFunction::ObjectiveFunction(MultSeqAlignment &alignment,
				     Phylogeny &phylogeny,
				     NthOrdRateMatrix &Q,
				     GradientMethod meth,
				     int dimensionality,
                                     GSL::Optimizer **optimizer,
                                     double gradientEpsilon,
                                     bool treeIsFixed,
                                     int numThreads,
				     MolecularSequenceType seqType,
				     ContextType contextType)
  : alignment(alignment), phylogeny(phylogeny), Q(Q),
    gradientMethod(meth), dimensionality(dimensionality),
    worstLikelihood(NEGATIVE_INFINITY), optimizer(optimizer),
    likelihood(NEGATIVE_INFINITY),
    dx(gradientEpsilon),
    treeIsFixed(treeIsFixed),
    numThreads(numThreads),
    seqType(seqType),
    contextType(contextType),
    numEvaluations(0),
    cacheHits(0)
{
  // ctor
}



bool ObjectiveFunction::pointIsValid(const GSL::Vector &point)
{
  int numMatrixParms=Q.numParameters()*Q.getNumMatrices();
  int numBranches=dimensionality-numMatrixParms;
  for(int i=0 ; i<numBranches ; ++i)
    if(!isFinite(point[i]) || point[i]<0) return false;
  for(int i=numBranches ; i<dimensionality ; ++i)
    {
      double x=point[i];
      if(!isFinite(x) || x<0 || x>1) return false;
    }
  return true;
}



double ObjectiveFunction::f(const GSL::Vector &currentPoint)
{
  if(!pointIsValid(currentPoint)) {
    //cout<<"leaving f(), returning NEG_INF"<<endl;
    return NEGATIVE_INFINITY;
  }
  if(cache.isDefined(currentPoint)) {
    ++cacheHits;
    return cache[currentPoint];
  }
  copyOut(currentPoint);
  NthOrdFelsenstein fel(contextType,phylogeny,&Q,alignment,
			seqType,numThreads);
  likelihood=fel.logLikelihood();
  ++numEvaluations;
  if(!isFinite(likelihood)) likelihood=worstLikelihood;
  else if(likelihood<worstLikelihood || !isFinite(worstLikelihood)) 
    worstLikelihood=likelihood;
  cout<<"     logL="<<likelihood<<" POINT="<<currentPoint<<endl;
  cache[currentPoint]=-likelihood;
  return -likelihood;
}



void ObjectiveFunction::gradient(const GSL::Vector &currentPoint,
				   GSL::Vector &gradient)
{
  if(!pointIsValid(currentPoint)) {
    //cout<<"setting gradient=NEGINF"<<endl;
    gradient.setAllTo(NEGATIVE_INFINITY);
    return;
  }
  gradient.resize(dimensionality);
  switch(gradientMethod)
    {
    case DIFFERENCING_METHOD:
      differencingMethod(currentPoint,gradient);
      break;
    case ANALYTICAL_METHOD:
      analyticalMethod(currentPoint,gradient);
      break;
    }
  //cout<<"gradient vector = "<<gradient<<endl;
  //cout<<Q<<endl;
  cout<<"\n#"<<((*optimizer)->iterationsUsed()+1)<<" logL="<<likelihood<<"\tGRAD="<<gradient.norm()<<endl;
}



void ObjectiveFunction::copyIn(GSL::Vector &v)
{
    int nextIndex=0;
    if(!treeIsFixed)
    {
        GetBranchLengths g(v);
        phylogeny.postorderTraversal(g);
        nextIndex=g.getNextIndex();
    }
    int n=Q.numParameters(), m=Q.getNumMatrices();
    for(int j=0 ; j<m ; ++j) {
      RateMatrix &M=Q.getIthMatrix(j);
      for(int i=0 ; i<n ; ++i, ++nextIndex)
        v[nextIndex]=M.getIthParm(i);
    }
}



void ObjectiveFunction::copyOut(const GSL::Vector &v)
{
    int nextIndex=0;
    if(!treeIsFixed)
    {
        SetBranchLengths s(v);
        phylogeny.postorderTraversal(s);
        nextIndex=s.getNextIndex();
    }
    int n=Q.numParameters(), m=Q.getNumMatrices();
    for(int j=0 ; j<m ; ++j) {
      RateMatrix &M=Q.getIthMatrix(j);
      for(int i=0 ; i<n ; ++i, ++nextIndex)
        M.setIthParm(i,v[nextIndex]);
      M.init();
    }
}



void ObjectiveFunction::differencingMethod(const GSL::Vector &currentPoint,
					   GSL::Vector &gradient)
{
  copyOut(currentPoint);
  NthOrdFelsenstein fel(contextType,phylogeny,&Q,alignment,seqType,numThreads);
  double L0=likelihood=fel.logLikelihood();
  //cout<<currentPoint<<endl;
  //cout<<"L0="<<L0<<endl;
  GSL::Vector perturbed=currentPoint;
  int n=gradient.getDim();
  //cout<<"n="<<n<<endl;
  for(int i=0 ; i<n ; ++i)
    {
        perturbed[i]+=dx;
        if(!pointIsValid(perturbed)) {
	  //cout<<i<<" not valid"<<endl;
	  perturbed[i]-=dx;
	  gradient[i]=0;
	  continue;
	}
        copyOut(perturbed);
	NthOrdFelsenstein fel(contextType,phylogeny,&Q,alignment,seqType,
			      numThreads);
        double L=fel.logLikelihood();
	//cout<<perturbed<<endl;
	//cout<<"L("<<i<<")="<<L<<endl;
        double dy=L-L0;
        double slope=dy/dx;
        gradient[i]=-slope;
        perturbed[i]-=dx;

	//###
	/*
        copyOut(currentPoint);
	NthOrdFelsenstein debug(contextType,phylogeny,&Q,alignment,seqType,
				numThreads);
        cout<<"DEBUG: "<<fel.logLikelihood()<<endl;
	GSL::Vector debugV;
	debugV.resize(currentPoint.getDim());
	copyIn(debugV);
	cout<<"DEBUGV: "<<debugV<<endl;
	*/
    }
}



void ObjectiveFunction::analyticalMethod(const GSL::Vector &currentPoint,
					   GSL::Vector &gradient)
{
  throw "Not implemented";
}



/****************************************************************
                    GetBranchLengths methods
 ****************************************************************/
GetBranchLengths::GetBranchLengths(GSL::Vector &v)
  : nextIndex(0), v(v)
{
  // ctor
}



void GetBranchLengths::processNode(InternalNode &node)
{
  v[nextIndex++]=node.getLeftDistance();
  v[nextIndex++]=node.getRightDistance();
}



void GetBranchLengths::processNode(RootNode &node)
{
  v[nextIndex++]=node.getBranchLength();
}



/****************************************************************
                    SetBranchLengths methods
 ****************************************************************/
SetBranchLengths::SetBranchLengths(const GSL::Vector &v)
  : nextIndex(0), v(v)
{
  // ctor
}



void SetBranchLengths::processNode(InternalNode &node)
{
  node.setLeftDistance(v[nextIndex++]);
  node.setRightDistance(v[nextIndex++]);
}



void SetBranchLengths::processNode(RootNode &node)
{
  node.setBranchLength(v[nextIndex++]);
}




