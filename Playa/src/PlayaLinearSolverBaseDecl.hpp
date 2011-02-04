/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_LINEARSOLVERBASEDECL_HPP
#define PLAYA_LINEARSOLVERBASEDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaSolverState.hpp"
#include "PlayaObjectWithVerbosity.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
  using namespace Teuchos;
  template <class Scalar>
  class LinearOperator;

  template <class Scalar>
  class Preconditioner;

  template <class Scalar>
  class PreconditionerFactory;

  template <class Scalar>
  class Vector;
  

  /** */
  template <class Scalar>
  class LinearSolverBase : public ObjectWithVerbosity
  {
  public:
    /** */
    LinearSolverBase(const ParameterList& params);

    /** */
    virtual ~LinearSolverBase(){;}

    /** */
    virtual SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
                                      const Vector<Scalar>& rhs,
                                      Vector<Scalar>& soln) const = 0;

    /** Change the convergence tolerance. Default does nothing. */
    virtual void updateTolerance(const double& tol) {;}

    /** Set a user-defined preconditioning operator. Default is an error. */
    virtual void setUserPrec(const PreconditionerFactory<Scalar>& pf);

    /** Set a user-defined preconditioning operator. Default is an error. */
    virtual void setUserPrec(const LinearOperator<Scalar>& P,
      const LinearSolver<Scalar>& pSolver);

    /** */
    const ParameterList& parameters() const ;

    /** */
    ParameterList& parameters();

    /** */
    std::string verbosityParam() const ;

    /** */
    template <typename T>
    static void setParameter(const ParameterList& params,
                             T* valuePtr, 
                             const std::string& paramName);
  private:
    ParameterList params_;
  };
}

#endif
