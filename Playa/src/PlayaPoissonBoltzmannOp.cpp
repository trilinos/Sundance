/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaPoissonBoltzmannOp.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Teuchos;

namespace Playa
{

PoissonBoltzmannOp::PoissonBoltzmannOp(int nLocal, const VectorType<double>& vecType)
  : NonlinearOperatorBase<double>(), J_(nLocal, vecType), importer_(),
    uLeftBC_(0.0), uRightBC_(2.0*log(cosh(1.0/sqrt(2.0))))
{
  setDomainAndRange(J_.domain(), J_.range());

  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();
  if (nProc > 1)
  {
    Array<int> ghosts;
    int low = J_.domain().baseGlobalNaturalIndex();
    int high = low + J_.domain().numLocalElements();
    if (rank != nProc - 1)
    {
      ghosts.append(high);
    }
    if (rank != 0) 
    {
      ghosts.append(low-1);
    }

    importer_ = vecType.createGhostImporter(J_.domain(), ghosts.size(), &(ghosts[0]));
  }
  else
  {
    importer_ = vecType.createGhostImporter(J_.domain(), 0, 0);
  }
}

Vector<double> PoissonBoltzmannOp::getInitialGuess() const
{
  Vector<double> rtn = J_.domain().createMember();

  rtn.setToConstant(0.5);

  return rtn;
}


LinearOperator<double> 
PoissonBoltzmannOp::computeJacobianAndFunction(Vector<double>& functionValue) const 
{
  Tabs tab;
  Out::os() << tab << "in PBOp::computeJacAndVec" << std::endl;
  J_.setEvalPoint(currentEvalPt());

  RCP<GhostView<double> > u;
  Out::os() << tab << "importing view" << std::endl;
  importer_->importView(currentEvalPt(), u);
  Out::os() << tab << "done importing view" << std::endl;
  int low = J_.domain().baseGlobalNaturalIndex();
  int high = low + J_.domain().numLocalElements();
  Out::os() << tab << "my indices are: " << low << ", " << high << std::endl;

  functionValue = J_.range().createMember();
  double h= J_.h();

  for (int r=low; r<high; r++)
  {
    Tabs tab1;
    double u_i = u->getElement(r);
    Out::os() << tab1 << r << " " << u_i << std::endl;
    double f = 0.0;
    if (r==0) 
    {
      f = u_i - uLeftBC_;
    }
    else if (r==J_.domain().dim()-1)
    {
      f = u_i - uRightBC_;
    }
    else
    {
      double u_plus = u->getElement(r+1);
      double u_minus = u->getElement(r-1);
      f = (u_plus + u_minus - 2.0*u_i)/h/h - exp(-u_i);
    }
    functionValue[r-low] = f;
  }

  Out::os() << tab << "done PBOp::computeJacAndVec" << std::endl;
  return J_.getOp();
}

Vector<double> PoissonBoltzmannOp::exactSoln() const
{
  Tabs tab;
  Out::os() << tab << "in PBOp::exactSoln" << std::endl;
  Vector<double> rtn = J_.domain().createMember();

  int low = J_.domain().baseGlobalNaturalIndex();
  int high = low + J_.domain().numLocalElements();

  double h= J_.h();
  
  for (int r=low; r<high; r++)
  {
    double x = r*h;
    double u = 2.0*log(cosh(x/sqrt(2.0)));
    rtn[r-low] = u;
  }

  Out::os() << tab << "done PBOp::exactSoln" << std::endl;
  return rtn;
}

}
