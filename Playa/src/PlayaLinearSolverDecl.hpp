/* @HEADER@ */
//
 /* @HEADER@ */


#ifndef PLAYA_LINEARSOLVERDECL_HPP
#define PLAYA_LINEARSOLVERDECL_HPP

#include "PlayaTabs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"
#include "Teuchos_TimeMonitor.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearSolverBaseImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

inline static Teuchos::Time& solveTimer() 
{
  static Teuchos::RCP<Teuchos::Time> rtn 
    = Teuchos::TimeMonitor::getNewTimer("linear solve"); 
  return *rtn;
}

namespace Playa
{
using namespace Teuchos;

  
/**
 *
 */
template <class Scalar>
class LinearSolver : public Playa::Handle<LinearSolverBase<Scalar> >
{
public:
  /** */
  LinearSolver() : Playa::Handle<LinearSolverBase<Scalar> >() {;}
  /** */
  LinearSolver( Playa::Handleable<LinearSolverBase<Scalar> >* rawPtr) 
    : Playa::Handle<LinearSolverBase<Scalar> >(rawPtr) {;}
  /** */
  LinearSolver(const RCP<LinearSolverBase<Scalar> >& smartPtr)
    : Playa::Handle<LinearSolverBase<Scalar> >(smartPtr) {;}


  /** Change the convergence tolerance. Default does nothing. */
  void updateTolerance(const double& tol) {this->ptr()->updateTolerance(tol);}

  /** Set a user-defined preconditioner */
  void setUserPrec(const LinearOperator<Scalar>& op,
    const LinearSolver<Scalar>& pSolver) ;

  /** Set a user-defined preconditioner */
  void setUserPrec(const PreconditionerFactory<Scalar>& pf);


  /** */
  SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
    const Vector<Scalar>& rhs,
    Vector<Scalar>& soln) const ;
    
    

  /** */
  const ParameterList& parameters() const ;

  /** */
  ParameterList& parameters() ;
};

  
template <class Scalar> inline 
SolverState<Scalar> LinearSolver<Scalar>
::solve(const LinearOperator<Scalar>& op,
  const Vector<Scalar>& rhs,
  Vector<Scalar>& soln) const
{
  Tabs tab;
  TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
    "null pointer in LinearSolver<Scalar>::solve()");

  TEST_FOR_EXCEPTION(rhs.ptr().get()==0, std::runtime_error,
    "null rhs pointer in LinearSolver<Scalar>::solve()");

  TEST_FOR_EXCEPTION(op.ptr().get()==0, std::runtime_error,
    "null op pointer in LinearSolver<Scalar>::solve()");

  TimeMonitor timer(solveTimer());

  PLAYA_MSG1(this->ptr()->verb(), 
    tab << "Solver(" << this->description() << ") starting solve");

  SolverState<Scalar> rtn = this->ptr()->solve(op, rhs, soln);

  PLAYA_MSG1(this->ptr()->verb(), 
    tab << "Solver(" << this->description() << ") done solve:");
  Tabs tab1;
  PLAYA_MSG2(this->ptr()->verb(), 
    tab << "state=" << rtn);

  return rtn;    
}

template <class Scalar> inline 
const ParameterList& LinearSolver<Scalar>::parameters() const 
{
  TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
    "null pointer in LinearSolver<Scalar>::parameters()");
  return this->ptr()->parameters();
}

template <class Scalar> inline 
ParameterList& LinearSolver<Scalar>::parameters() 
{
  TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
    "null pointer in LinearSolver<Scalar>::parameters()");
  return this->ptr()->parameters();
}

  

  

}

#endif
