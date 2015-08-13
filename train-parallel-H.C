/****************************************************************
 train-parallel-H.C (H=higher order)
 bmajoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3.
 ****************************************************************/
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/MultiAlignment.H"
#include "BOOM/Vector.H"
#include "BOOM/Random.H"
#include "BOOM/Constants.H"
#include "BOOM/ConfigFile.H"
#include "BOOM/GSL/Optimizer.H"
#include "BOOM/Environment.H"
#include "BOOM/Time.H"
#include "BOOM/Map.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/Exceptions.H"
#include "BOOM/BitSet.H"
#include "BOOM/AlphabetMap.H"
#include "BOOM/PureDnaAlphabet.H"
#include "Phylogeny.H"
#include "NthOrdRateMatrix.H"
#include "RateMatrixType.H"
#include "NthOrdFelsenstein.H"
#include "NmerFelsenstein.H"
#include "Message.H"
#include "AlignmentNmerTable.H"
using namespace std;
using namespace BOOM;


const double pi=2.0*acos(0.0);

inline double mapToProbability(double x)
{
  return atan(x)/pi+0.5;
}

inline double mapFromProbability(double x) 
{
  return tan(pi*(x-0.5));
}

GSL::Vector mapToProbabilities(const GSL::Vector &v,int numBranches)
{
  GSL::Vector r=v;
  int dimensionality=r.getDim();

  // Constrain branch lengths to be >0
  for(int i=0 ; i<numBranches ; ++i)
    if(r[i]<0.01) {r[i]=0.01;}

  // Constrain matrix parameters to be between 0 and 1
  for(int i=numBranches ; i<dimensionality ; ++i)
    r[i]=mapToProbability(r[i]);

  // Return mapped point
  return r;
}

GSL::Vector mapFromProbabilities(const GSL::Vector &v,int numBranches)
{
  GSL::Vector r=v;
  int dimensionality=r.getDim();
  for(int i=numBranches ; i<dimensionality ; ++i)
    r[i]=mapFromProbability(r[i]);
  return r;
}



/****************************************************************
                       enum GradientMethod
 ****************************************************************/
enum GradientMethod {
  DIFFERENCING_METHOD,
  ANALYTICAL_METHOD
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
  int numProgressiveReps;
  bool fixOneBranchLength;
  GSL::Optimizer *optimizer;
  String optimizerTypeStr, stoppingCriterionStr;
  double tolerance, gradientThreshold, gradientEpsilon, stepSize;
  bool treeIsFixed;
  int numProcesses, processID;
  int order, maxIterations;
  bool dual, progressive;
  Alphabet *alphabet;
  BitSet gapSymbols;
  MolecularSequenceType seqType;
  float noiseSigma; // for generating random (higher-order) rate matrix
  ContextType contextType;
  GSL::Vector optimalPoint, lowerBounds, upperBounds;
  double optimalValue; // log-likelihood of optimalPoint
  Array1D<double> equilibriumFreqs;
  GradientMethod gradientMethod;
  bool regularize;
  double regularizationParm;
  int maxSeqLen; // maximum (total) seq length for training
  AlignmentNmerTable *nmerTable;
  AlphabetMap *alphabetMap; // maps gapped alphabet to ungapped alphabet
  int numBranches;

  int master(int argc,char *argv[]);
  int slave(int argc,char *argv[]);
  void getEquilibriumFreqs(MultSeqAlignment &,
			   Array1D<double> &equilibriumFreqs);
  void optimize();
  void initBounds();
  void progressiveZerothOrder();
  void copyIn(GSL::Vector &v) {copyIn(v,*Q,*phylogeny,treeIsFixed);}
  void copyOut(GSL::Vector &v) {copyOut(v,*Q,*phylogeny,treeIsFixed);}
public:
  static GSL::Vector originalBranchLengths;
  static double branchCoefficient; // scaling factor for branch lengths

  Application();
  int main(int argc,char *argv[]);
  static void copyIn(GSL::Vector &,NthOrdRateMatrix &,Phylogeny &,
		     bool treeIsFixed);
  static void copyOut(const GSL::Vector &,NthOrdRateMatrix &,Phylogeny &,
		      bool treeIsFixed);
};


GSL::Vector Application::originalBranchLengths; // decl of static member
double Application::branchCoefficient=1.0;



/****************************************************************
                     class ObjectiveFunction
 ****************************************************************/
class ObjectiveFunction : public GSL::ObjectiveFunction
{
  bool fixOneBranchLength;
  GradientMethod gradientMethod;
  MultSeqAlignment &alignment;
  Phylogeny &phylogeny;
  NthOrdRateMatrix &Q;
  int dimensionality, numProcesses, numBranches;
  double worstLikelihood, likelihood;
  GSL::Optimizer **optimizer;
  double dx;
  bool treeIsFixed;
  MolecularSequenceType seqType;
  ContextType contextType;
  BOOM::Map<const GSL::Vector,double> cache;
  int cacheHits;
  int numEvaluations;
  double firstScore, bestScore;
  bool regularize;
  double regularizationParm;
  GSL::Vector bestPoint;
  
  void differencingMethod(const GSL::Vector &currentPoint,
			  GSL::Vector &gradient);
  void twoPointAsymmetric(const GSL::Vector &currentPoint,
			  GSL::Vector &gradient);
  void twoPointSymmetric(const GSL::Vector &currentPoint,
			  GSL::Vector &gradient);
  void fourPointSymmetric(const GSL::Vector &currentPoint,
			  GSL::Vector &gradient);
  void analyticalMethod(const GSL::Vector &currentPoint,
			GSL::Vector &gradient);
  bool pointIsValid(const GSL::Vector &point);
  void truncate(GSL::Vector &point);
  GSL::Vector constrain(const GSL::Vector &);
  double constrained_f(const GSL::Vector &currentPoint);
public:
  ObjectiveFunction(MultSeqAlignment &,Phylogeny &,NthOrdRateMatrix &,
		    GradientMethod,int dimensionality,GSL::Optimizer **,
		    double gradientEpsilon,bool treeIsFixed,int numProcesses,
		    bool fixOneBranchLength,
		    MolecularSequenceType=DNA,ContextType=CT_RCO,
		    bool regularize=false,double regularizationParm=1);
  virtual double f(const GSL::Vector &currentPoint);
  virtual void gradient(const GSL::Vector &currentPoint,
			GSL::Vector &gradient);
  int getNumEvaluations() const {return numEvaluations;}
  int getNumCacheHits() const {return cacheHits;}
  void copyIn(GSL::Vector &v) 
    {Application::copyIn(v,Q,phylogeny,treeIsFixed);}
  void copyOut(const GSL::Vector &v) 
    {Application::copyOut(v,Q,phylogeny,treeIsFixed);}
  double getFirstScore() {return firstScore;}
  double getBestScore() {return bestScore;}
  GSL::Vector &getBestPoint() {return bestPoint;}
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
    catch(const String &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(const RootException &e)
      {
	cerr << "exception: "<< e.getMessage() << endl;
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
      treeIsFixed(false), noiseSigma(0.05), 
      equilibriumFreqs(DnaDashDotAlphabet::global().size()),
      gapSymbols(256), nmerTable(NULL)
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD,&processID);
  int exitCode;

  // Determine whether this process is a MASTER or a SLAVE
  if(processID==0) {
    exitCode=master(argc,argv);
    MPI_Abort(MPI_COMM_WORLD,exitCode);
    MPI_Finalize();
  }
  else {
    exitCode=slave(argc,argv);
    MPI_Finalize();
  }

  // Clean up & exit
  //MPI_Finalize();
  //MPI_Abort(MPI_COMM_WORLD,0);
  return exitCode;
}



/****************************************************************
                             MASTER
 ****************************************************************/
int Application::master(int argc,char *argv[])
{
    Time stopwatch;
    stopwatch.startCounting();

    // Process command line
    CommandLine cmd(argc,argv,"m:tdpM:");
    if(cmd.numArgs()!=5)
      throw String(
"\ntrain-parallel-H <in.phy> <in.maf> <out.phy> <out.matrix> <config-file>\n\
\n\
options: \n\
         -m matrix = seed the matrix rather than initializing randomly\n\
         -M matrix = seed matrix, but optimize at next higher order only\n\
         -t = don't modify the tree parameters\n\
         -d = use dual contexts\n\
         -p = (progressive) optimize lower orders before higher ones\n\
\n"
);
    String phyFile=cmd.arg(0);
    String mafFile=cmd.arg(1);
    String outPhy=cmd.arg(2);
    String outMatrix=cmd.arg(3);
    String configFilename=cmd.arg(4);
    treeIsFixed=cmd.option('t');
    randomize();
    dual=cmd.option('d');
    progressive=cmd.option('p');

    cout.precision(8);

    // Process the configuration file
    cout<<"MASTER processing config file"<<endl;
    ConfigFile configFile(configFilename);
    String Mtype=configFile.lookupOrDie("matrix-type");
    maxIterations=configFile.getIntOrDie("max-iterations");
    optimizerTypeStr=configFile.lookupOrDie("optimizer");
    stoppingCriterionStr=configFile.lookupOrDie("stopping-criterion");
    tolerance=configFile.getDoubleOrDie("tolerance");
    gradientThreshold=configFile.getDoubleOrDie("gradient-threshold");
    gradientEpsilon=configFile.getDoubleOrDie("gradient-epsilon");
    stepSize=configFile.getDoubleOrDie("step-size");
    int maxOrder=configFile.getIntOrDie("order");
    seqType=seqTypeFromString(configFile.lookupOrDie("sequence-type"));
    if(configFile.isDefined("noise-sigma")) 
      noiseSigma=configFile.getFloatOrDie("noise-sigma");
    contextType=contextTypeFromString(configFile.lookupOrDie("context-type"));
    order=progressive ? 0 : maxOrder;
    gradientMethod=DIFFERENCING_METHOD;
    regularize=configFile.getBoolOrDie("regularize");
    regularizationParm=
      (regularize ? configFile.getDoubleOrDie("regularization-parm") : 1);
    maxSeqLen=configFile.getIntOrDie("max-seq-length");
    numProgressiveReps=configFile.getIntOrDie("num-progressive-replicates");
    fixOneBranchLength=configFile.getBoolOrDie("fix-one-branch-length");

    // Load the alignments
    cout<<"MASTER loading alignments"<<endl;
    Vector<MultiAlignment*> alignments;
    MultiAlignment::loadMAF(mafFile,alignments);
    MultiAlignment *alignment=MultiAlignment::combine(alignments,true);//###
    cout<<"\tMASTER converting to uppercase...\n";
    alignment->toupper();
    cout<<"\tMASTER acquiring gap symbols...\n";
    alphabet=&DnaDashDotAlphabet::global();
    gapSymbols.addMember(alphabet->lookup('-'));
    gapSymbols.addMember(alphabet->lookup('.'));
    gapSymbols.addMember(alphabet->lookup('N'));
    cout<<"\tMASTER converting to symbolic alignment...\n";
    MultSeqAlignment &A=*new MultSeqAlignment(*alignment,*alphabet,gapSymbols);
    this->alignment=&A;
    cout<<"\MASTER tacquiring alphabet map...\n";
    alphabetMap=new DropGapMapping(*alphabet,PureDnaAlphabet::global());
    //### should delete alignment now

    // Load the phylogeny
    cout<<"MASTER loading phylogeny"<<endl;
    phylogeny=new Phylogeny(phyFile);
    phylogeny->attachAlignment(A);
    numBranches=phylogeny->getNumNodes()-1;
    originalBranchLengths.resize(numBranches);
    GetBranchLengths g(originalBranchLengths);
    phylogeny->postorderTraversal(g);
    cout<<*phylogeny<<endl; // ### DEBUGGING

    // Delete alignment columns having a gap in the root sequence
    cout<<"MASTER deleting target gaps"<<endl;
    int rootID=phylogeny->getRoot()->getID();
    A.deleteTargetGaps(rootID);
    if(contextType==CT_NMER) {
      cout<<"constructing nmer table"<<endl;
      nmerTable=new AlignmentNmerTable(A,*alphabetMap,order);
    }

    // Get equilibrium frequencies
    cout<<"MASTER getting EQ freqs"<<endl;
    getEquilibriumFreqs(A,equilibriumFreqs);
    cout<<equilibriumFreqs<<endl; // ### DEBUGGING

    // Load or generate the rate matrix
    cout<<"MASTER loading or generating matrix"<<endl;
    if(cmd.option('m')) 
      {
	Q=NthOrdRateMatrix::load(cmd.optParm('m'));
	matrixType=Q->getMatrixType();
	order=Q->getOrder();
      }
    else if(cmd.option('M')) 
      {
	Q=NthOrdRateMatrix::load(cmd.optParm('M'));
	matrixType=Q->getMatrixType();
	order=Q->getOrder();
      }
    else 
      {
	int ord=dual ? 2*order : order;
	matrixType=rateMatrixTypeFromString(Mtype);
	Q=NthOrdRateMatrix::random(dual,noiseSigma,seqType,ord,
				   matrixType,equilibriumFreqs);
      }
    int numMatrices=Q->getNumMatrices();
    int numParms=Q->numParameters();
    dimensionality=numMatrices*numParms;
    dimensionality+=(treeIsFixed ? 1 : numBranches);

    // Initialize the slaves
    cout<<"MASTER init slaves"<<endl;
    int numCols=A.getLength();
    if(numCols>maxSeqLen) numCols=maxSeqLen;
    float colsPerSlave=numCols/(numProcesses-1.0);
    int begin=0;
    float end=colsPerSlave;
    Message message(dimensionality);
    void *buffer=message.getBuffer();
    int bufferSize=message.getBufferSize();
    const int TAG=0;
    for(int slaveID=1 ; slaveID<numProcesses ; ++slaveID) {
      if(slaveID+1>=numProcesses) end=float(numCols);
      int iEnd=int(end);
      message.packInterval(begin,iEnd);
      MPI_Send(buffer,bufferSize,MPI_BYTE,slaveID,TAG,MPI_COMM_WORLD);
      begin=iEnd;
      end+=colsPerSlave;
    }

    // Optimize
    cout<<"MASTER optimizer"<<endl;
    optimalPoint.resize(dimensionality);
    copyIn(optimalPoint);
    optimalValue=NEGATIVE_INFINITY;
    if(cmd.option('M')) progressive=true;
    else if(progressive && order==0) progressiveZerothOrder();
    else optimize();

    // option "-t" : allow branch lengths to vary; re-optimize
    if(treeIsFixed) {
      message.packTreeNotFixed();
      for(int slaveID=1 ; slaveID<numProcesses ; ++slaveID)
	MPI_Send(buffer,bufferSize,MPI_BYTE,slaveID,TAG,MPI_COMM_WORLD);
      treeIsFixed=false;
      dimensionality=Q->getNumMatrices()*Q->numParameters()+numBranches;
      message.resize(dimensionality);
      buffer=message.getBuffer();
      bufferSize=message.getBufferSize();
      optimize();
    }

    // option "-p" : incrementally raise the order & re-optimize
    if(progressive) {
      int finalOrder=dual ? 2*maxOrder : maxOrder;
      int orderIncrement=dual ? 2 : 1;
      int initOrder=order+1;//###(cmd.option('m') ? order+1 : orderIncrement)
      for(order=initOrder ; order<=finalOrder ; order+=orderIncrement) {
	if(contextType==CT_NMER)
	  nmerTable=nmerTable->nextHigherOrder();
	message.packNextHigherOrder(order,optimalPoint);
	for(int slaveID=1 ; slaveID<numProcesses ; ++slaveID)
	  MPI_Send(buffer,bufferSize,MPI_BYTE,slaveID,TAG,MPI_COMM_WORLD);
	NthOrdRateMatrix *newQ=Q->nextHigherOrder();
	newQ->averageLowerOrderModel();
	delete Q;
	Q=newQ;
	numMatrices=Q->getNumMatrices();
	numParms=Q->numParameters();
	dimensionality=numMatrices*numParms;
	dimensionality+=(treeIsFixed ? 1 : numBranches);
	message.resize(dimensionality);
	buffer=message.getBuffer();
	bufferSize=message.getBufferSize();
	/*
	{
	  // DEBUGGING:
	  ObjectiveFunction f(A,*phylogeny,*Q,gradientMethod,
			      dimensionality,&optimizer,gradientEpsilon,
			      treeIsFixed,numProcesses,fixOneBranchLength,
			      seqType,contextType,
			      regularize,regularizationParm);
	  GSL::Vector pt;
	  pt.resize(dimensionality);
	  f.copyIn(pt);
	  double x=f.f(pt);
	}
	*/
	optimize();
      }
    }

    // Generate output
    ofstream os(outPhy.c_str());
    phylogeny->printOn(os);
    Q->save(outMatrix);

    // Terminate the slaves
    message.packQuit();
    for(int slaveID=1 ; slaveID<numProcesses ; ++slaveID) {
      MPI_Send(buffer,bufferSize,MPI_BYTE,slaveID,TAG,MPI_COMM_WORLD);
    }

    stopwatch.stopCounting();
    cout<<"MASTER Done.  Elapsed time: "<<stopwatch.elapsedSeconds()<<" sec ("
	<<numProcesses<<" CPUs)"<<endl;
    return 0;
  }


/*
  This function performs several optimizations at zeroth order to find the
  best zeroth order initial model.  Each replicate starts over with a random
  point.
 */
void Application::progressiveZerothOrder() {
  // Get the initial point
  GSL::Vector initialPoint;
  initialPoint.resize(dimensionality);
  copyIn(initialPoint);

  // Iteratively try out different random points for the zeroth-order model
  int numReplicates=numProgressiveReps;
  GSL::Vector bestPoint;
  double bestLL;
  for(int i=0 ; i<numReplicates ; ++i) {
    // Run the optimizer
    optimize();

    // See if we have a new winner
    if(i==0 || optimalValue>bestLL) {
      bestPoint=optimalPoint;
      bestLL=optimalValue;
      cout<<"NEW OPTIMUM: "<<bestLL<<endl;
    }
    
    // Get ready for next iteration (try another random point)
    if(i+1<numReplicates) {
      copyOut(initialPoint);
      int ord=dual ? 2*order : order;
      delete Q;
      Q=NthOrdRateMatrix::random(dual,noiseSigma,seqType,ord,
				 matrixType,equilibriumFreqs);
    }
  }
  
  optimalPoint=bestPoint;
  optimalValue=bestLL;
  cout<<"ZEROTH ORDER OPTIMUM: "<<bestLL<<" @ "<<bestPoint<<endl;
  copyOut(optimalPoint);
}



/****************************************************************
                              SLAVE
 ****************************************************************/
int Application::slave(int argc,char *argv[])
{
    // Process command line
    CommandLine cmd(argc,argv,"m:tdpM:");
    if(cmd.numArgs()!=5)
      throw String(
"\ntrain-parallel-H <in.phy> <in.maf> <out.phy> <out.matrix> <config-file>\n\
\n\
options: \n\
         -m matrix = seed the matrix rather than initializing randomly\n\
         -t = don't modify the tree parameters\n\
         -d = use dual contexts\n\
\n"
);
    String phyFile=cmd.arg(0);
    String mafFile=cmd.arg(1);
    String configFilename=cmd.arg(4);
    treeIsFixed=cmd.option('t');
    dual=cmd.option('d');
    progressive=cmd.option('p');
    
    // Process the configuration file
    //cout<<"SLAVE: processing config file" <<endl;
    ConfigFile configFile(configFilename);
    String Mtype=configFile.lookupOrDie("matrix-type");
    maxIterations=configFile.getIntOrDie("max-iterations");
    optimizerTypeStr=configFile.lookupOrDie("optimizer");
    stoppingCriterionStr=configFile.lookupOrDie("stopping-criterion");
    tolerance=configFile.getDoubleOrDie("tolerance");
    gradientThreshold=configFile.getDoubleOrDie("gradient-threshold");
    gradientEpsilon=configFile.getDoubleOrDie("gradient-epsilon");
    stepSize=configFile.getDoubleOrDie("step-size");
    seqType=seqTypeFromString(configFile.lookupOrDie("sequence-type"));
    int maxOrder=configFile.getIntOrDie("order");
    if(configFile.isDefined("noise-sigma")) 
      noiseSigma=configFile.getFloatOrDie("noise-sigma");
    contextType=contextTypeFromString(configFile.lookupOrDie("context-type"));
    fixOneBranchLength=configFile.getBoolOrDie("fix-one-branch-length");
    order=progressive ? 0 : maxOrder;
    gradientMethod=DIFFERENCING_METHOD;

    // Load the alignments
    //cout<<"SLAVE: loading alignments"<<endl;
    Vector<MultiAlignment*> alignments;
    MultiAlignment::loadMAF(mafFile,alignments);
    MultiAlignment *alignment=MultiAlignment::combine(alignments,true);
    //cout<<"\tSLAVE: converting to uppercase"<<endl;
    alignment->toupper();
    alphabet=&DnaDashDotAlphabet::global();
    //cout<<"\tSLAVE: acquiring gap symbols"<<endl;
    gapSymbols.addMember(alphabet->lookup('-'));
    gapSymbols.addMember(alphabet->lookup('.'));
    gapSymbols.addMember(alphabet->lookup('N'));
    //cout<<"\tSLAVE: converting to symbolic alignment"<<endl;
    MultSeqAlignment &A=*new MultSeqAlignment(*alignment,*alphabet,gapSymbols);
    this->alignment=&A;
    alphabetMap=new DropGapMapping(*alphabet,PureDnaAlphabet::global());

    // Load the phylogeny
    //cout<<"SLAVE: loading phylogeny"<<endl;
    phylogeny=new Phylogeny(phyFile);
    phylogeny->attachAlignment(A);
    numBranches=phylogeny->getNumNodes()-1;
    originalBranchLengths.resize(numBranches);
    GetBranchLengths g(originalBranchLengths);
    phylogeny->postorderTraversal(g);

    // Delete alignment columns having a gap in the root sequence
    //cout<<"SLAVE: deleting target gaps"<<endl;
    int rootID=phylogeny->getRoot()->getID();
    A.deleteTargetGaps(rootID);
    if(contextType==CT_NMER) {
      //cout<<"SLAVE: constructing nmer table"<<endl;
      nmerTable=new AlignmentNmerTable(A,*alphabetMap,order);
    }

    // Get equilibrium frequencies
    //cout<<"SLAVE: getting EQ frequencies"<<endl;
    getEquilibriumFreqs(A,equilibriumFreqs);

    // Load or generate the rate matrix
    //cout<<"SLAVE: loading or generating rate matrix"<<endl;
    if(cmd.option('m')) 
      {
	Q=NthOrdRateMatrix::load(cmd.optParm('m'));
	matrixType=Q->getMatrixType();
      }
    else if(cmd.option('M'))
      {
	Q=NthOrdRateMatrix::load(cmd.optParm('M'));
	matrixType=Q->getMatrixType();
      }
    else 
      {
	int ord=dual ? 2*order : order;
	matrixType=rateMatrixTypeFromString(Mtype);
	Q=NthOrdRateMatrix::random(dual,noiseSigma,seqType,ord,
				   matrixType,equilibriumFreqs);
      }
    dimensionality=Q->getNumMatrices()*Q->numParameters();
    dimensionality+=(treeIsFixed ? 1 : numBranches);

    // Serve the master...
    //cout<<"SLAVE: waiting for instructions from master..."<<endl;
    int intervalBegin, intervalEnd;
    Message message(dimensionality);
    void *buffer=message.getBuffer();
    int bufferSize=message.getBufferSize();
    int sourceID=0;
    MPI_Status status;
    double L;
    bool done=false;
    const TAG=0;
    const MASTER=0;
    while(!done) {
      MPI_Recv(buffer,bufferSize,MPI_BYTE,sourceID,MPI_ANY_TAG,
	       MPI_COMM_WORLD,&status);
      switch(message.unpack()) 
	{
	case MSG_SET_INTERVAL:
	  //cout<<"SLAVE SET INTERVAL"<<endl;
	  intervalBegin=message.getBegin();
	  intervalEnd=message.getEnd();
	  break;
	case MSG_EVALUATE:
	  {
	    //cout<<"SLAVE EVALUATE"<<endl;
	    copyOut(message.getParms());
	    //cout<<"SLAVE copyout"<<endl;
	    Q->averageLowerOrderModel();
	    //cout<<"SLAVE avelowerorder"<<endl;
	    ContextType ct=contextType;
	    if(ct==CT_NMER) {
	      NmerFelsenstein F(*phylogeny,Q,A,nmerTable);
	      L=F.logLikelihood(intervalBegin,intervalEnd);
	    }
	    else {
	      //cout<<"SLAVE allocate fels"<<endl;
	      NthOrdFelsenstein F(ct,*phylogeny,Q,A,seqType,nmerTable);
	      //cout<<"SLAVE compute likelihood"<<endl;
	      L=F.logLikelihood(intervalBegin,intervalEnd);
	    }
	    //cout<<"SLAVE packlikelihood"<<endl;
	    message.packLikelihood(L);
	    //cout<<"SLAVE mpi-send"<<endl;
	    MPI_Send(buffer,bufferSize,MPI_BYTE,MASTER,TAG,MPI_COMM_WORLD);
	    //cout<<"SLAVE: sent"<<endl;
	  }
	  break;
	case MSG_TREE_NOT_FIXED:
	  treeIsFixed=false;
	  dimensionality=Q->getNumMatrices()*Q->numParameters()+numBranches;
	  message.resize(dimensionality);
	  buffer=message.getBuffer();
	  bufferSize=message.getBufferSize();
	  break;
	case MSG_NEXT_HIGHER_ORDER: 
	  {
	    if(contextType==CT_NMER)
	      nmerTable=nmerTable->nextHigherOrder();
	    order=message.getOrder();
	    copyOut(message.getParms());
	    NthOrdRateMatrix *newQ=Q->nextHigherOrder();
	    delete Q;
	    Q=newQ;
	    dimensionality=Q->getNumMatrices()*Q->numParameters();
	    dimensionality+=treeIsFixed ? 1 : numBranches;
	    message.resize(dimensionality);
	    buffer=message.getBuffer();
	    bufferSize=message.getBufferSize();
	  }
	  break;
	case MSG_QUIT:
	  //cout<<"SLAVE QUIT"<<endl;
	  done=true;
	  break;
	default:
	  cout<<"SLAVE: unknown message!"<<endl;
	  break;
	}
    }

    return 0;
  }



void Application::getEquilibriumFreqs(MultSeqAlignment &alignment,
				      Array1D<double> &equilibriumFreqs)
{
  int n=alignment.getNumTracks();
  equilibriumFreqs.resize(alphabet->size()); // ACGTN-.
  equilibriumFreqs.setAllTo(0.0);
  for(int i=0 ; i<n ; ++i)
    {
      AlignmentSeq &track=alignment.getIthTrack(i);
      const Sequence &seq=track.getSeq();
      int L=seq.getLength();
      for(int i=0 ; i<L ; ++i) {
	Symbol a=seq[i];
	char c=alphabet->lookup(a);
	switch(c) {
	case 'A':
	case 'C':
	case 'G':
	case 'T':
	  ++equilibriumFreqs[a]; 
	}
      }
    }
  int sampleSize=0;
  int numAlpha=alphabet->size();
  for(Symbol s=0 ; s<numAlpha ; ++s) {
    char c=alphabet->lookup(s);
    switch(c) {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
      sampleSize+=equilibriumFreqs[s];
    }
  }
  for(Symbol s=0 ; s<numAlpha ; ++s)
    equilibriumFreqs[s]/=sampleSize;
}



void Application::initBounds() {
  lowerBounds.resize(dimensionality);
  upperBounds.resize(dimensionality);
  lowerBounds.setAllTo(0.0001);
  upperBounds.setAllTo(0.99);

  // The first N parameters are branch lengths, the rest are matrix parms:
  int numMatrixParms=Q->numParameters()*Q->getNumMatrices();
  int numBranches=dimensionality-numMatrixParms;
  for(int i=0 ; i<numBranches ; ++i) upperBounds[i]=1000;

  if(fixOneBranchLength) {
    // Require one of the branch lengths to have length 1.0
    upperBounds[0]=1.0;
    lowerBounds[0]=1.0;
  }
}


void Application::optimize()
{
  ObjectiveFunction f(*alignment,*phylogeny,*Q,gradientMethod,
		      dimensionality,&optimizer,gradientEpsilon,
                      treeIsFixed,numProcesses,fixOneBranchLength,
		      seqType,contextType,
		      regularize,regularizationParm);
  GSL::Vector initialPoint;
  initialPoint.resize(dimensionality);
  f.copyIn(initialPoint);
  cout<<"INITIAL POINT: "<<initialPoint<<endl;
  GSL::OptimizerType optimizerType=
      GSL::stringToOptimizerType(optimizerTypeStr);
  GSL::StoppingCriterion stoppingCriterion=
      GSL::stringToStoppingCriterion(stoppingCriterionStr);

  cout<<"init bounds"<<endl;
  initBounds();
  
  cout<<"mapping initial point"<<endl;
  initialPoint=mapFromProbabilities(initialPoint,numBranches);
  cout<<"allocating optimizer"<<endl;
  optimizer=new GSL::Optimizer(optimizerType,f,initialPoint,stepSize,
                               stoppingCriterion,tolerance,gradientThreshold,
                               maxIterations);
  /*optimizer=new GSL::ConstrainedOptimizer(optimizerType,f,initialPoint,
					  lowerBounds,
					  upperBounds,stepSize,
					  stoppingCriterion,
					  likelihoodThreshold,maxIterations);*/
  cout<<"running optimizer"<<endl;
  optimizer->run();
  cout<<"optimization complete."<<endl;

  // Collect the optimal point
  GSL::Vector bestPoint;
  bestPoint=f.getBestPoint();
  double bestValue=f.getBestScore();
  if(!isFinite(optimalValue) || bestValue>optimalValue 
     || bestPoint.getDim()!=optimalPoint.getDim()) // ### 7/23/07
    {
      optimalValue=bestValue;
      optimalPoint=bestPoint;
    }
  delete optimizer;
  copyOut(optimalPoint);
  cout<<"BEST POINT for this iteration: "<<bestPoint<<endl;
  cout<<"logL increase: "
      <<f.getFirstScore()
      <<" to "
      <<bestValue
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
                                     int numProcesses,
				     bool fixOneBranchLength,
				     MolecularSequenceType seqType,
				     ContextType contextType,
				     bool regularize,
				     double regularizationParm)
  : alignment(alignment), phylogeny(phylogeny), Q(Q),
    gradientMethod(meth), dimensionality(dimensionality),
    worstLikelihood(NEGATIVE_INFINITY), optimizer(optimizer),
    likelihood(NEGATIVE_INFINITY),
    dx(fabs(gradientEpsilon)),
    treeIsFixed(treeIsFixed),
    numProcesses(numProcesses),
    seqType(seqType),
    contextType(contextType),
    numEvaluations(0),
    cacheHits(0),
    firstScore(NEGATIVE_INFINITY),
    bestScore(NEGATIVE_INFINITY),
    regularize(regularize),
    regularizationParm(regularizationParm),
    fixOneBranchLength(fixOneBranchLength)
{
  int numMatrixParms=Q.numParameters()*Q.getNumMatrices();
  numBranches=dimensionality-numMatrixParms;
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



void ObjectiveFunction::truncate(GSL::Vector &point)
{
    int corrections=0;
    for(int i=0 ; i<dimensionality ; ++i) {
      double &x=point[i];
      if(x<=0) {x=dx; ++corrections;}
      else if(x>=1) {x=1-dx; ++corrections;}
    }
    if(corrections>0) cout<<"TRUNCATED IN "<<corrections<<" DIMENSIONS\n";
}



GSL::Vector ObjectiveFunction::constrain(const GSL::Vector &v)
{
  return mapToProbabilities(v,numBranches);
}



double ObjectiveFunction::f(const GSL::Vector &v)
{
  return constrained_f(constrain(v));
}



double ObjectiveFunction::constrained_f(const GSL::Vector &currentPt)
{
  // PRECONDITION: currentPt has already been constrained

  GSL::Vector currentPoint=currentPt;
  if(cache.isDefined(currentPoint)) {
    ++cacheHits;
    return cache[currentPoint];
  }

  // Send work orders to the slaves
  Message message(dimensionality);
  void *buffer=message.getBuffer();
  int bufferSize=message.getBufferSize();
  MPI_Status status;
  message.packEvaluate(currentPoint);
  const int TAG=0;
  for(int i=1 ; i<numProcesses ; ++i)
    MPI_Send(buffer,bufferSize,MPI_BYTE,i,TAG,MPI_COMM_WORLD);
  
  // Wait for the slaves to return their likelihoods
  likelihood=0;
  for(int i=1 ; i<numProcesses ; ++i) {
    MPI_Recv(buffer,bufferSize,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,
	     MPI_COMM_WORLD,&status);
    MessageType msg=message.unpack();
    likelihood+=message.getLikelihood();
  }
  ++numEvaluations;
  
  // Perform error checking, and return the likelihood
  if(!isFinite(likelihood)) likelihood=worstLikelihood;
  else if(likelihood<worstLikelihood || !isFinite(worstLikelihood)) 
    worstLikelihood=likelihood;
  cout<<"    logL="<<likelihood<<" POINT="<<currentPoint<<endl;

  if(!isFinite(firstScore)) {
    firstScore=bestScore=likelihood;
    bestPoint=currentPt;
  }
  else if(likelihood>bestScore || currentPt.getDim()>bestPoint.getDim()) {
    bestScore=likelihood;
    bestPoint=currentPt;
  }

  double answer=-likelihood;
  if(regularize) {
    double norm=currentPt.norm();
    double reg2=regularizationParm*regularizationParm;
    double penalty=norm*norm/(2*reg2);
    answer+=penalty;
    cout<<"\t\tregularization term: "<<penalty<<endl;
  }
  cache[currentPoint]=answer;
  return answer;
}



void ObjectiveFunction::gradient(const GSL::Vector &currentPt,
				   GSL::Vector &gradient)
{
  // PRECONDITION: currentPt has *not* been constrained

  static iterations=0;
  if(gradient.getDim()!=dimensionality) throw "gradient has wrong dimension";
  //gradient.resize(dimensionality); // ### 12/24/07
  
  switch(gradientMethod)
    {
    case DIFFERENCING_METHOD:
      differencingMethod(currentPt,gradient);
      break;
    case ANALYTICAL_METHOD:
      analyticalMethod(currentPt,gradient);
      break;
    }
  cout<<"#"<<(++iterations)<<" logL="<<likelihood<<"\tGRAD="<<gradient.norm()<<endl;
}



void ObjectiveFunction::differencingMethod(const GSL::Vector &currentPt,
					   GSL::Vector &gradient)
{
  // PRECONDITION: currentPt has *not* been constrained

  //twoPointAsymmetric(currentPt,gradient);
  twoPointSymmetric(currentPt,gradient);
  //fourPointSymmetric(currentPt,gradient);
}


void ObjectiveFunction::twoPointAsymmetric(const GSL::Vector &currentPt,
					   GSL::Vector &gradient)
{
  // Requires only D+1 evaluations, for dimensionality D

  // PRECONDITION: currentPt has *not* been constrained

  GSL::Vector currentPoint=currentPt;
  double L0=f(currentPoint);
  GSL::Vector perturbed=currentPoint;
  int n=gradient.getDim();
  int firstIndex;
  if(fixOneBranchLength) {
    gradient[0]=0;
    firstIndex=1;
  }
  else firstIndex=0;
  for(int i=firstIndex ; i<n ; ++i)
    {
      perturbed[i]+=dx;
      double L=f(perturbed);
      double dy=L-L0;
      double slope=dy/dx;
      gradient[i]=slope;
      perturbed[i]-=dx;
    }
}



void ObjectiveFunction::twoPointSymmetric(const GSL::Vector &currentPt,
					  GSL::Vector &gradient)
{
  // Requires 2D evaluations, for dimensionality D

  // PRECONDITION: currentPt has *not* been constrained

  GSL::Vector currentPoint=currentPt;
  GSL::Vector perturbed=currentPoint;
  int n=gradient.getDim();
  int firstIndex;
  if(fixOneBranchLength) {
    gradient[0]=0;
    firstIndex=1;
  }
  else firstIndex=0;
  const double twoDX=2*dx;
  for(int i=firstIndex ; i<n ; ++i)
    {
      double &pi=perturbed[i];
      pi-=dx;
      double L0=f(perturbed);
      pi+=twoDX;
      double L=f(perturbed);
      double dy=L-L0;
      double slope=dy/twoDX;
      gradient[i]=slope;
      pi-=dx;
    }
}



void ObjectiveFunction::fourPointSymmetric(const GSL::Vector &currentPt,
					   GSL::Vector &gradient)
{
  // Requires only D+1 evaluations, for dimensionality D

  // PRECONDITION: currentPt has *not* been constrained

  GSL::Vector currentPoint=currentPt;
  GSL::Vector perturbed=currentPoint;
  int n=gradient.getDim();
  int firstIndex;
  if(fixOneBranchLength) {
    gradient[0]=0;
    firstIndex=1;
  }
  else firstIndex=0;
  const double twoDX=2*dx, twelveDX=12*dx;
  for(int i=firstIndex ; i<n ; ++i)
    {
      double &pi=perturbed[i];
      pi-=twoDX;
      double sum=f(perturbed);
      pi+=dx;
      sum-=8*f(perturbed);
      pi+=twoDX;
      sum+=8*f(perturbed);
      pi+=dx;
      sum-=f(perturbed);
      gradient[i]=sum/twelveDX;
      pi-=twoDX;
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




/****************************************************************
                      global functions
 ****************************************************************/

void Application::copyIn(GSL::Vector &v,NthOrdRateMatrix &Q,
			 Phylogeny &phylogeny,bool treeIsFixed)
{
  int nextIndex;
  if(treeIsFixed) {
    v[0]=branchCoefficient;
    nextIndex=1;
  }
  else {
    GetBranchLengths g(v);
    phylogeny.postorderTraversal(g);
    nextIndex=g.getNextIndex();
  }

  int n=Q.numParameters(), m=Q.getNumMatrices();
  if(Q.getType()==MT_NMER) {
    NmerRateMatrix &M=static_cast<NmerRateMatrix&>(Q);
    for(int i=0 ; i<n ; ++i, ++nextIndex)
      v[nextIndex]=M.getIthParm(i);    
  }
  else {
    for(int j=0 ; j<m ; ++j) {
      RateMatrix &M=Q.getIthMatrix(j);
      for(int i=0 ; i<n ; ++i, ++nextIndex)
	v[nextIndex]=M.getIthParm(i);
    }
  }
}



void Application::copyOut(const GSL::Vector &v,NthOrdRateMatrix &Q,
			  Phylogeny &phylogeny,bool treeIsFixed)
{
  int nextIndex;
  if(treeIsFixed) {
    branchCoefficient=v[0];
    GSL::Vector lengths=originalBranchLengths;
    lengths.scale(branchCoefficient);
    SetBranchLengths s(lengths);
    phylogeny.postorderTraversal(s);
    nextIndex=1;
  }
  else {
    SetBranchLengths s(v);
    phylogeny.postorderTraversal(s);
    nextIndex=s.getNextIndex();
  }

  int n=Q.numParameters(), m=Q.getNumMatrices();
  if(Q.getType()==MT_NMER) {
    NmerRateMatrix &M=static_cast<NmerRateMatrix&>(Q);
    for(int i=0 ; i<n ; ++i, ++nextIndex) {
      M.setIthParm(i,v[nextIndex]);
    }
  }
  else {
    for(int j=0 ; j<m ; ++j) {
      RateMatrix &M=Q.getIthMatrix(j);
      for(int i=0 ; i<n ; ++i, ++nextIndex)
	M.setIthParm(i,v[nextIndex]);
      M.init();
    }
  }
}



