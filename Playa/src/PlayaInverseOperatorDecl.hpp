/* @HEADER@ */
//
 /* @HEADER@ */

#ifndef PLAYA_INVERSEOPERATOR_DECL_HPP
#define PLAYA_INVERSEOPERATOR_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "Teuchos_RCP.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaSolverState.hpp"

namespace Playa
{
using Teuchos::RCP;

/** 
 * PlayaInverseOperator represents the inverse of some other operator.  An
 * inverse operator object will contain an operator and a solver.  The 
 * operator data member is the operator whose inverse this represents.  The
 * solver data member is the solver that will be used in applying the
 * inverse.  If the solver is null, the operator is assumed to have
 * self-contained ability to solve systems, as for a dense matrix that 
 * does solves by factoring and backsolves.
 */
template <class Scalar> 
class InverseOperator : public LinearOpWithSpaces<Scalar>,
                        public Printable
{
public:
  /**
   * Ctor with a linear operator and a solver specified.
   */
  InverseOperator(const LinearOperator<Scalar>& op, 
    const LinearSolver<Scalar>& solver);

  /** Virtual dtor */
  virtual ~InverseOperator(){;}

  /** 
   * Apply the operator. 
   * 
   * \param applyType Indicates whether to apply the operator, its transpose,
   * or its conjugate transpose. 
   * \param in The vector on which the operator is to act
   * \param out The vector into which the result of the operation 
   * is to be written. This vector should already be initialized by the
   * appropriate space.
   **/
  virtual void apply(
    Teuchos::ETransp applyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const ;


  
  /** */
  void print(std::ostream& os) const ;

  /** */
  LinearOperator<Scalar> op() const {return op_;}


private:
  const LinearOperator<Scalar> op_;
  const LinearSolver<Scalar> solver_;  
  std::string msg_;
};


/** \brief Implicit inverse operator. */
template <class Scalar> 
LinearOperator<Scalar> 
inverse(const LinearOperator<Scalar>& op, 
  const LinearSolver<Scalar>& solver);
  

}

#endif
