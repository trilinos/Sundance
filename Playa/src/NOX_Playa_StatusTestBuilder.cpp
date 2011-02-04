// $Id$ 
// $Source$ 

//@HEADER
//   
//@HEADER

#include "NOX_Playa_StatusTestBuilder.hpp"         
#include "NOX_StatusTest_NormF.H"         
#include "NOX_StatusTest_NormUpdate.H"         
#include "NOX_StatusTest_SafeCombo.hpp"         
#include "NOX_StatusTest_MaxIters.H"         
#include "Teuchos_TestForException.hpp"   

using namespace NOX;
using namespace NOX::NOXPlaya;
using namespace Teuchos;

RCP<StatusTest::Generic> 
StatusTestBuilder::makeStatusTest(const ParameterList& params)
{
  TEST_FOR_EXCEPTION(!params.isSublist("Status Test"), runtime_error,
                     "did not find Status Test sublist in " << params);

  ParameterList testSublist = params.sublist("Status Test");

  double fTol = 1.0e-15;
  double dxTol = 1.0e-15;
  int maxiters = 20;
  if (testSublist.isParameter("Tolerance"))
    {
      fTol = getParameter<double>(testSublist, "Tolerance");
    }
  if (testSublist.isParameter("Residual Tolerance"))
    {
      fTol = getParameter<double>(testSublist, "Residual Tolerance");
    }
  if (testSublist.isParameter("Step Tolerance"))
    {
      dxTol = getParameter<double>(testSublist, "Step Tolerance");
    }
  if (testSublist.isParameter("Max Iterations"))
    {
      maxiters = getParameter<int>(testSublist, "Max Iterations");
    }

  RCP<StatusTest::Generic> A = rcp(new StatusTest::NormF(fTol));
  RCP<StatusTest::Generic> B = rcp(new StatusTest::MaxIters(maxiters));
  RCP<StatusTest::Generic> C = rcp(new StatusTest::NormUpdate(dxTol));
  RCP<StatusTest::Generic> AB 
    = rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, A, B));
  RCP<StatusTest::Generic> ABC 
    = rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, AB, C));
  
  return ABC;
}



