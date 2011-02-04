//@HEADER@

//@HEADER@

#include "BelosBlockGmresSolMgr.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "PlayaBelosAdapter.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaPreconditioner.hpp"
#include "PlayaPreconditionerFactory.hpp"
#include "PlayaParameterListPreconditionerFactory.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaSimpleDiagonalOpImpl.hpp"
#include "PlayaSimpleIdentityOpImpl.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaAztecSolver.hpp"
#include "PlayaAnasaziEigensolverDecl.hpp"
#include "PlayaAnasaziAdapter.hpp"
#include "PlayaEigensolver.hpp"
#include "PlayaMatrixLaplacian1D.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaInverseOperatorDecl.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "AnasaziMVOPTester.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaInverseOperatorImpl.hpp"
#include "PlayaAnasaziEigensolverImpl.hpp"

#endif

using namespace Teuchos;
using namespace Playa;
using namespace PlayaOps;


int main(int argc, char *argv[]) 
{
  typedef Teuchos::ScalarTraits<double> ST;

  try
  {
    GlobalMPISession session(&argc, &argv);

    MPIComm::world().synchronize();

    VectorType<double> type = new EpetraVectorType();


    ParameterXMLFileReader reader("anasazi-ml.xml");
    ParameterList solverParams = reader.getParameters().sublist("Eigensolver");

    /* create the range space  */
    int nLocalRows = 40;
    MatrixLaplacian1D builder(nLocalRows, type);
    typedef Anasazi::MultiVec<double> MV;
    typedef Anasazi::Operator<double> OP;

    LinearOperator<double> A = builder.getOp();
    LinearOperator<double> M;

    Teuchos::RCP<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::BasicOutputManager<double>() );
    MyOM->setVerbosity(Anasazi::Warnings);

    int nv = 1;
    RCP<const Array<Vector<double> > > initMV = rcp(new Array<Vector<double> >(nv));
    RCP<Array<Vector<double> > > nc 
      = rcp_const_cast<Array<Vector<double> > >(initMV);
    for (int i=0; i<nv; i++) 
    {
      (*nc)[i] = A.domain().createMember();
      (*nc)[i].randomize();
    }

#ifdef BLARF
    bool mvPass = Anasazi::TestMultiVecTraits<double,Array<Vector<double> > >(MyOM,initMV);
    if (mvPass) Out::os() << "******* MV unit test PASSED ******* " << endl;
    else Out::os() << "******* MV unit test FAILED ******* " << endl;

    RCP<const LinearOperator<double> > APtr = rcp(&A, false);
    bool opPass = Anasazi::TestOperatorTraits<double, Array<Vector<double> >, LinearOperator<double>  >(MyOM, initMV, APtr);
    if (opPass) Out::os() << "******* OP unit test PASSED ******* " << endl;
    else Out::os() << "******* OP unit test FAILED ******* " << endl;
#endif
    Eigensolver<double> solver = new AnasaziEigensolver<double>(solverParams);
    
    Array<Vector<double> > ev;
    Array<std::complex<double> > ew;
    
    solver.solve(A, M, ev, ew);

    Out::os() << "Eigenvalues are " << ew << endl;

    const double pi = 4.0*atan(1.0);
    double err = 0.0;
    for (int i=0; i<ev.size(); i++)
    {
      double x = (i+1)*pi;
      err += ::fabs(ew[i].real()-x*x)/x/x;
    }
    err = err / ew.size();
    
    Out::os() << "error = " << err << endl;
    if (err < 0.01)
    {
      cout << "Belos poisson solve test PASSED" << std::endl;
      return 0;
    }
    else
    {
      cout << "Belos poisson solve test FAILED" << std::endl;
      return 1;
    }
  }
  catch(std::exception& e)
  {
    cout << "Caught exception: " << e.what() << std::endl;
    return -1;
  }
}

