/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUnaryMinusEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceUnaryMinus.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;




UnaryMinusEvaluator
::UnaryMinusEvaluator(const UnaryMinus* expr,
                      const EvalContext& context)
  : UnaryEvaluator<UnaryMinus>(expr, context)
{
  int vecResultIndex = 0;
  int constResultIndex = 0;
  
  for (int i=0; i<sparsity()->numDerivs(); i++)
    {
      /* Determine the index into which the result will be written */
      bool resultIsConstant = sparsity()->state(i)==ConstantDeriv; 

      if (!resultIsConstant)
        {
          addVectorIndex(i, vecResultIndex);
          vecResultIndex++;
        }
      else
        {
          addConstantIndex(i, constResultIndex);
          constResultIndex++;
        }
    }
}

void UnaryMinusEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RefCountPtr<EvalVector> >& vectorResults) const
{
  TimeMonitor timer(evalTimer());
  Tabs tab;
  SUNDANCE_OUT(verbosity() > VerbLow,
               tab << "------- UnaryMinusEvaluator -------");


  /* evaluate the argument */
  Array<RefCountPtr<EvalVector> > argVectorResults;
  Array<double> argConstantResults;

  evalOperand(mgr, constantResults, vectorResults);


  if (verbosity() > VerbLow)
    {
      cerr << tab << "operand results" << endl;
      argSparsity()->print(cerr, vectorResults,
                           constantResults);
    }

  for (unsigned int i=0; i<constantResults.size(); i++)
    {
      constantResults[i] *= -1;
    }

  for (unsigned int i=0; i<vectorResults.size(); i++)
    {
      vectorResults[i]->multiply_S(-1.0);
    }

  
  
}


