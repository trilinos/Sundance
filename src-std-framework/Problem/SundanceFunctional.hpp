/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FUNCTIONAL_H
#define SUNDANCE_FUNCTIONAL_H

#include "SundanceDefs.hpp"
#include "SundanceLinearProblem.hpp"
#include "SundanceNonlinearProblem.hpp"
#include "TSFNonlinearOperator.hpp"
#include "TSFLinearSolver.hpp"
#include "TSFVector.hpp"
#include "TSFVectorType.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Teuchos;

  /**
   *
   */
  class Functional
  {
  public:
    /** */
    Functional(const Mesh& mesh, const Expr& integral, 
               const TSFExtended::VectorType<double>& vecType);

    /** */
    Functional(const Mesh& mesh, const Expr& integral, 
               const Expr& essentialBC,
               const TSFExtended::VectorType<double>& vecType);

    /** */
    LinearProblem linearVariationalProb(const Expr& var,
                                        const Expr& varEvalPts,
                                        const Expr& unk,
                                        const Expr& fixed,
                                        const Expr& fixedEvalPts) const ;

    
    /** */
    NonlinearOperator<double>
    nonlinearVariationalProb(const Expr& var,
                             const Expr& varEvalPts,
                             const Expr& unk,
                             const Expr& unkEvalPts,
                             const Expr& fixed,
                             const Expr& fixedEvalPts) const ;


    double evaluate() const ;

  private:
    Mesh mesh_;

    Expr integral_;

    Expr bc_;

    TSFExtended::VectorType<double> vecType_;
    
  };
}


#endif
