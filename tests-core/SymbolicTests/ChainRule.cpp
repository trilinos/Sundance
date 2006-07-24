#include "SundanceChainRuleEvaluator.hpp"
#include "SundanceUnknownFunctionStub.hpp"

using namespace SundanceUtils;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;
using SundanceCore::List;

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

			Expr u = new UnknownFunctionStub("u");
			Expr v = new UnknownFunctionStub("v");
			Expr w = new UnknownFunctionStub("w");


      MultiSet<int> Q0 = makeMultiSet(0);
      MultiSet<int> Q00 = makeMultiSet(0, 0);
      MultiSet<int> Q000 = makeMultiSet(0, 0, 0);
      MultiSet<int> Q01 = makeMultiSet(0, 1);
      MultiSet<int> Q012 = makeMultiSet(0, 1, 2);
      MultiSet<int> Q001 = makeMultiSet(0, 0, 1);

      MultipleDeriv uv = makeDeriv(u,v);
      MultipleDeriv uvw = makeDeriv(u,v,w);
      MultipleDeriv uvv = makeDeriv(u,v,v);
      MultipleDeriv uu = makeDeriv(u,u);
      MultipleDeriv uuu = makeDeriv(u,u,u);

      cout << "uu/0 = " << ChainRuleEvaluator::chainRuleBins(uu, Q0) << endl;
      cout << "uv/0 = " << ChainRuleEvaluator::chainRuleBins(uv, Q0) << endl;
      cout << "uvv/0 = " << ChainRuleEvaluator::chainRuleBins(uvv, Q0) << endl;
      cout << "uvw/0 = " << ChainRuleEvaluator::chainRuleBins(uvw, Q0) << endl;
      cout << "uuu/0 = " << ChainRuleEvaluator::chainRuleBins(uuu, Q0) << endl;

      cout << "uu/00 = " << ChainRuleEvaluator::chainRuleBins(uu, Q00) << endl;
      cout << "uv/00 = " << ChainRuleEvaluator::chainRuleBins(uv, Q00) << endl;
      cout << "uvv/00 = " << ChainRuleEvaluator::chainRuleBins(uvv, Q00) << endl;
      cout << "uvw/00 = " << ChainRuleEvaluator::chainRuleBins(uvw, Q00) << endl;

      cout << "uu/01 = " << ChainRuleEvaluator::chainRuleBins(uu, Q01) << endl;
      cout << "uv/01 = " << ChainRuleEvaluator::chainRuleBins(uv, Q01) << endl;
      cout << "uvv/01 = " << ChainRuleEvaluator::chainRuleBins(uvv, Q01) << endl;
      cout << "uvw/01 = " << ChainRuleEvaluator::chainRuleBins(uvw, Q01) << endl;

      cout << "uuu/000 = " 
           << ChainRuleEvaluator::chainRuleBins(uuu, Q000) << endl;
      cout << "uuu/001 = " 
           << ChainRuleEvaluator::chainRuleBins(uuu, Q001) << endl;
      cout << "uuu/012 = " 
           << ChainRuleEvaluator::chainRuleBins(uuu, Q012) << endl;

      cout << "uvv/000 = " 
           << ChainRuleEvaluator::chainRuleBins(uvv, Q000) << endl;
      cout << "uvv/001 = " 
           << ChainRuleEvaluator::chainRuleBins(uvv, Q001) << endl;
      cout << "uvv/012 = " 
           << ChainRuleEvaluator::chainRuleBins(uvv, Q012) << endl;

      cout << "uvw/000 = " 
           << ChainRuleEvaluator::chainRuleBins(uvw, Q000) << endl;
      cout << "uvw/001 = " 
           << ChainRuleEvaluator::chainRuleBins(uvw, Q001) << endl;
      cout << "uvw/012 = " 
           << ChainRuleEvaluator::chainRuleBins(uvw, Q012) << endl;

      
    }
	catch(exception& e)
		{
			Out::println(e.what());
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
