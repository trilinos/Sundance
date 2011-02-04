/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_BLOCKTRIANGULARSOLVER_DECL_HPP
#define PLAYA_BLOCKTRIANGULARSOLVER_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverDecl.hpp" 


namespace Playa
{
/** */
template <class Scalar>
class BlockTriangularSolver : public LinearSolverBase<Scalar>,
                              public Playa::Handleable<LinearSolverBase<Scalar> >
{
public:
  /** */
  BlockTriangularSolver(const LinearSolver<Scalar>& solver);

  /** */
  BlockTriangularSolver(const Array<LinearSolver<Scalar> >& solvers);

  /** */
  virtual ~BlockTriangularSolver(){}

  /** */
  virtual SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
    const Vector<Scalar>& rhs,
    Vector<Scalar>& soln) const ;

  /* */
  GET_RCP(LinearSolverBase<Scalar>);
private:
  Array<LinearSolver<Scalar> > solvers_;
};

}

#endif
