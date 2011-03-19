/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaSteepestDescent.hpp"
#include "PlayaOptBuilder.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaOut.hpp"
#include "PlayaLinearCombinationImpl.hpp"

#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"

#include <fstream>

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaVectorImpl.hpp"
#include "PlayaBlockIteratorImpl.hpp"
#endif

using namespace Playa;
using std::ofstream;

#ifdef BLARF
class Rosenbrock : public ObjectiveBase
{
public:
  Rosenbrock(int n, double alpha, const VectorType<double>& vecType);

  void  evalGrad(const Vector<double>& x, double& f, 
    Vector<double>& grad) const ;

  void eval(const Vector<double>& x, double& f) const ;

  Vector<double> getInit() const ;

private:
  int n_;
  double alpha_;
};

void Rosenbrock::eval(const Vector<double>& x, double& f) const
{
  f = 0.0;
  for (int i=0; i<n_; i+=2)
  {
    double p = x[i+1] - x[i]*x[i];
    double q = 1.0-x[i];
    f += alpha_ * p*p + q*q;
  }
}
#endif


class Ellipsoid : public ObjectiveBase
{
public:
  Ellipsoid(int n, const VectorType<double>& vecType)
    : n_(n), vs_(vecType.createEvenlyPartitionedSpace(MPIComm::self(), n))
    {}

  void  evalGrad(const Vector<double>& x, double& f, 
    Vector<double>& grad) const ;

  void eval(const Vector<double>& x, double& f) const ;

  Vector<double> getInit() const ;

private:
  int n_;
  VectorSpace<double> vs_;
};

void Ellipsoid::eval(const Vector<double>& x, double& f) const
{
  f = 0.0;
  for (int i=0; i<n_; i++)
  {
    f += ::sqrt(2+i)*x[i]*x[i];
  }
}

void Ellipsoid::evalGrad(const Vector<double>& x, double& f, 
  Vector<double>& grad) const
{
  f = 0.0;
  for (int i=0; i<n_; i++)
  {
    f += ::sqrt(i+2)*x[i]*x[i];
    grad[i] = 2.0*::sqrt(i+2)*x[i];
  }
}

Vector<double> Ellipsoid::getInit() const
{
  Vector<double> rtn = vs_.createMember();
  rtn.setToConstant(1.0);
  rtn[0] = 0.5;
  return rtn;
}


int main(int argc, char *argv[])
{
  int stat = 0;
  try
  {
    GlobalMPISession session(&argc, &argv);

    VectorType<double> vecType = new EpetraVectorType();
    double testTol = 1.0e-4;

    
    int n = 200;
    RCP<ObjectiveBase> ell = rcp(new Ellipsoid(n, vecType));

    Array<RCP<ObjectiveBase> > probs = tuple(ell);
    Array<string> algs 
      = tuple<string>("basicLMBFGS", "steepestDescent");

    int numFail = 0;
    int testCount = 0;

    for (int i=0; i<probs.size(); i++)
    {
      RCP<ObjectiveBase> obj = probs[i];
      Vector<double> xInit = obj->getInit();
      obj->fdCheck(xInit, 0);

      for (int j=0; j<algs.size(); j++)
      {
        RCP<UnconstrainedOptimizerBase> opt 
          = OptBuilder::createOptimizer(algs[j] + ".xml");

        RCP<ConvergenceMonitor> mon = rcp(new ConvergenceMonitor());

        opt->setVerb(0);

        OptState state = opt->run(obj, xInit, mon);

        if (state.status() != Opt_Converged)
        {
          Out::root() << "optimization failed: " << state.status() << endl;
          numFail++;
        }
        else
        {
          Out::root() << "optimization succeeded!" << endl ;
          
          Vector<double> exactAns = xInit.copy();
          exactAns.zero();
          
          double locErr = (exactAns - state.xCur()).norm2();
          Out::root() << "error in location = " << locErr << endl;
          if (locErr > testTol) 
          {
            Out::root() << "test FAILED" << endl;
            numFail++;
          }
          
          string monName = algs[j] + ".dat";
          ofstream os(monName.c_str());
          mon->write(os);
        }
      }
      testCount++;
    }
    
    if (numFail > 0)
    {
      Out::root() << "detected " << numFail 
                  << " FAILURES out of " << testCount << " tests" << endl;
      stat = -1;
    }
    else
    {
      Out::root() << "all tests PASSED" << endl;
    }
  }
  catch(std::exception& e)
  {
    std::cerr << "Caught exception: " << e.what() << endl;
    stat = -1;
  }
  return stat;
}
