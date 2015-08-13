/****************************************************************
 eigen-decompose.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/GSL/Matrix.H"
#include "BOOM/GSL/Vector.H"
#include "RateMatrix.H"
#include "BOOM/Random.H"
using namespace std;
using namespace BOOM;


class Application
{
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
    CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=1)
      throw String("eigen-decompose <matrix-file>");
    String infile=cmd.arg(0);

    RateMatrix M(infile,DNA);

    for(int trial=0 ; trial<9000000 ; ++trial)
      {
	// Generate the equilibrium frequencies
	double pi[4], sum=0;
	for(int i=0 ; i<4 ; ++i)
	  sum+=(pi[i]=Random0to1());
	for(int i=0 ; i<4 ; ++i)
	  pi[i]/=sum;
	
	// Generate transition/transvertion rates
	double alpha=Random0to1();
	double beta=1-alpha;
	double kappa=Random0to1();
	double chi=Random0to1();
	double omega=Random0to1();
	double tau=Random0to1();
	
	// Populate the rate matrix
	double a[4][4]={
	  0,           beta*pi[1],  alpha*pi[2], chi*pi[3],
	  beta*pi[0],  0,           kappa*pi[2], omega*pi[3],
	  alpha*pi[0], kappa*pi[1], 0,           tau*pi[3],
	  chi*pi[0],   omega*pi[1], tau*pi[2],   0
	};
	
	// Adjust the diagonal entries to maintain a sum of zero
	for(int i=0 ; i<4 ; ++i)
	  {
	    double sum=a[i][0]+a[i][1]+a[i][2]+a[i][3];
	    a[i][i]=-sum;
	  }
	
	GSL::Matrix m(4,4);
	for(int i=0 ; i<4 ; ++i)
	  for(int j=0 ; j<4 ; ++j)
	    m(i,j)=a[j][i];
	
	GSL::Matrix eigenvectors;
	GSL::Vector eigenvalues;
	GSL::Vector imagParts;
	m.getEigenVectors(eigenvectors,eigenvalues,imagParts);
	
	for(int i=0 ; i<4 ; ++i)
	  if(fabs(imagParts[i])>1e-15)
	    cout<<m<<"\n I["<<i<<"]="<<imagParts[i];
      }

    /*
    GSL::Matrix e1=m;
    for(int i=0 ; i<1000 ; ++i)
      {
	GSL::Matrix temp;
	e1.times(m,temp);
	e1=temp;
	GSL::Vector freq;
	e1.getRow(0,freq);
	double P=
	  freq[0]+freq[1]+freq[2]+freq[3];
	freq.scale(1/P);
      }
    cout<<e1<<endl;
    */

    /*
    for(int i=0 ; i<4 ; ++i)
      {
	GSL::Vector eigenvector;
	eigenvectors.getColumn(i,eigenvector);
	double P=
	  eigenvector[0]+eigenvector[1]+eigenvector[2]+eigenvector[3];
	//eigenvector.scale(1/P);
	cout<<eigenvector<<endl;
	double lambda=eigenvalues[i];
	GSL::Vector Ax, lambdaX;
	m.multiply(eigenvector,Ax);
	lambdaX=eigenvector;
	lambdaX.scale(lambda);
	cout<<lambda<<endl<<Ax<<endl<<lambdaX<<endl<<endl;
      }
    */

    /*
    cout<<"---------------------------"<<endl;
    for(int i=0 ; i<4 ; ++i)
      for(int j=0 ; j<4 ; ++j)
	{
	  GSL::Vector vi, vj;
	  eigenvectors.getColumn(i,vi);
	  eigenvectors.getColumn(j,vj);
	  cout<<i<<","<<j<<" => "<<vi.dotProduct(vj)<<endl;

	}
    */

    return 0;
  }

