/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFunctional.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceAssembler.hpp"

using namespace SundanceStdFwk;
using namespace SundanceCore;
using namespace SundanceStdMesh;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;


Functional::Functional(const Mesh& mesh, const Expr& integral, 
                       const TSFExtended::VectorType<double>& vecType)
  : mesh_(mesh),
    integral_(integral),
    bc_(),
    vecType_(vecType)
{
  
}

Functional::Functional(const Mesh& mesh, const Expr& integral,  
                       const Expr& essentialBC,
                       const TSFExtended::VectorType<double>& vecType)
  : mesh_(mesh),
    integral_(integral),
    bc_(essentialBC),
    vecType_(vecType)
{
  
}


LinearProblem Functional::linearVariationalProb(const Expr& var,
                                                const Expr& varEvalPts,
                                                const Expr& unk,
                                                const Expr& fixed,
                                                const Expr& fixedEvalPts) const
{

  Array<Expr> zero(unk.size());
  for (int i=0; i<unk.size(); i++) 
    {
      Expr z = new ZeroExpr();
      zero[i] = z;
    }

  Expr unkEvalPts = new ListExpr(zero);

  RefCountPtr<EquationSet> eqn 
    = rcp(new EquationSet(integral_, bc_, var, varEvalPts,
                          unk, unkEvalPts, fixed, fixedEvalPts));

  RefCountPtr<Assembler> assembler 
    = rcp(new Assembler(mesh_, eqn, vecType_));

  return LinearProblem(assembler);
}

NonlinearOperator<double> Functional
::nonlinearVariationalProb(const Expr& var,
                           const Expr& varEvalPts,
                           const Expr& unk,
                           const Expr& unkEvalPts,
                           const Expr& fixed,
                           const Expr& fixedEvalPts) const
{
  RefCountPtr<EquationSet> eqn 
    = rcp(new EquationSet(integral_, bc_, var, varEvalPts,
                          unk, unkEvalPts, fixed, fixedEvalPts));

  RefCountPtr<Assembler> assembler 
    = rcp(new Assembler(mesh_, eqn, vecType_));

  return new NonlinearProblem(assembler, unkEvalPts);
}

FunctionalEvaluator Functional::evaluator(const Expr& var,
                                          const Expr& varEvalPts,
                                          const Expr& fixed,
                                          const Expr& fixedEvalPts) const 
{
  return FunctionalEvaluator(mesh_, integral_, bc_,
                             var, 
                             varEvalPts, 
                             fixed, 
                             fixedEvalPts,
                             vecType_);
}


FunctionalEvaluator Functional::evaluator(const Expr& var,
                                          const Expr& varEvalPts) const 
{
  return FunctionalEvaluator(mesh_, integral_, bc_,
                             var, 
                             varEvalPts,
                             vecType_);
}

