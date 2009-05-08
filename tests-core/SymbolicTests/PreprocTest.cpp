#include "SundanceSymbPreprocessor.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "Teuchos_TestingHelpers.hpp"

using namespace SundanceUtils;
using namespace SundanceCore;
using namespace Teuchos;
using namespace TSFExtended;
using SundanceCore::List;


#define TEST_THROW(code, passFail) \
  TEUCHOS_TEST_THROW( code, std::exception, Out::os(), passFail)

#define TEST_NOTHROW(code, passFail) \
  TEUCHOS_TEST_NOTHROW( code, Out::os(), passFail)

bool validateFuncTypeChecking()
{
  Expr ux = new UnknownFunctionStub("ux");
  Expr vx = new TestFunctionStub("vx");
  Expr uy = new UnknownFunctionStub("uy");
  Expr vy = new TestFunctionStub("vy");
  Expr uz = new UnknownFunctionStub("uz");
  Expr vz = new TestFunctionStub("vz");

  Expr v = List(vx,vy,vz);
  Expr u = List(ux,uy,uz);

  Expr mixup = List(vx, uy, vz); // mix of test & unknown
  Expr dup = List(vx, vx, vz); // list with duplicates

  bool passFail = true;
  /* */
  Out::os() << "Testing detection of mixed-up function types" << endl;
  TEST_THROW(
    SymbPreprocessor::processInputFuncs<UnknownFuncElement>(mixup, makeZeros(v)), passFail);
  
  /* */
  Out::os() << "Testing detection of duplicated functions" << endl;
  TEST_THROW(
    SymbPreprocessor::processInputFuncs<SymbolicFuncElement>(dup, makeZeros(v)), passFail);
  
  /* */
  Out::os() << "Testing detection of invalid evaluation points" << endl;
  TEST_THROW(
    SymbPreprocessor::processInputFuncs<SymbolicFuncElement>(v, u), passFail);
  

  /* */
  Out::os() << "Testing processing of good input" << endl;
  TEST_NOTHROW(
    SymbPreprocessor::processInputFuncs<SymbolicFuncElement>(v, makeZeros(v)), passFail);
  
  return passFail;
} 

int main(int argc, char** argv)
{
  bool pass = true;
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      pass = pass && validateFuncTypeChecking();
    }
	catch(exception& e)
		{
      pass = false;
			Out::println(e.what());
		}

  if (pass)
  {
    Out::os() << "test PASSED" << endl;
    return 0;
  }
  else 
  {
    Out::os() << "test FAILED" << endl;
    return -1;
  }

  
}
