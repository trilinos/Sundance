/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_LINEARPROBLEM_H
#define SUNDANCE_LINEARPROBLEM_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "TSFLinearOperator.hpp"
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
  class LinearProblem 
    : public TSFExtended::ObjectWithVerbosity<LinearProblem>
  {
  public:
    /** Empty ctor */
    LinearProblem();
    
    /** Construct with a mesh, equation set, bcs, test and unknown funcs,
     * and a vector type */
    LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
                  const Expr& test, const Expr& unk, 
                  const TSFExtended::VectorType<double>& vecType);
    
    /** */
    Expr solve(const LinearSolver<double>& solver) const ;

    /** */
    Vector<double> getRHS() const ;

    /** */
    LinearOperator<double> getOperator() const ;


  private:
      
    /** */
    RefCountPtr<Assembler> assembler_;

    /** */
    mutable LinearOperator<double> A_;

    /** */
    mutable Vector<double> rhs_;
    
  };
}


#endif
