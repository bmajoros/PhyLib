
#include <iostream>
#include "BOOM/GSL/GammaDistribution.H"
#include "BOOM/CommandLine.H"
using namespace std;
using namespace BOOM;
using namespace GSL;

const double DELTA=0.01;
const double MAX=1000;

/****************************************************************
                       class Application
 ****************************************************************/
class Application
{
  CommandLine *cmd;
public:
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
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }




int Application::main(int argc,char *argv[])
{
  // Process command line
  cmd=new CommandLine(argc,argv,"");
  if(cmd->numArgs()!=2) throw String("\ngamma <shape> <scale>\n\n");
  double shape=cmd->arg(0);
  double scale=cmd->arg(1);

  GammaDistribution gamma(shape,scale);
  for(double i=DELTA ; i<MAX ; i+=DELTA) {
    double p=gamma.probabilityEquals(i);
    if(p<0.0001) break;
    cout<<i<<"\t"<<p<<endl;
  }
}
