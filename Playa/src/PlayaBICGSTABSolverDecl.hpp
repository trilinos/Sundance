/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_BICGSTABSOLVER_DECL_HPP
#define PLAYA_BICGSTABSOLVER_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaKrylovSolver.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "PlayaOut.hpp"
#include "Teuchos_Describable.hpp"


namespace Playa
{
using namespace Teuchos;

/**
 *
 */
template <class Scalar>
class BICGSTABSolver : public KrylovSolver<Scalar>,
                       public Playa::Handleable<LinearSolverBase<Scalar> >,
                       public Printable,
                       public Describable
{
public:
  /** */
  BICGSTABSolver(const ParameterList& params = ParameterList());

  /** */
  BICGSTABSolver(const ParameterList& params,
    const PreconditionerFactory<Scalar>& precond);

  /** */
  virtual ~BICGSTABSolver(){;}

  /** \name Printable interface */
  //@{
  /** Write to a stream  */
  void print(std::ostream& os) const ;
  //@}
    
  /** \name Describable interface */
  //@{
  /** Write a brief description */
  std::string description() const {return "BICGSTABSolver";}
  //@}

  /** \name Handleable interface */
  //@{
  /** Return a ref count pointer to a newly created object */
  virtual RCP<LinearSolverBase<Scalar> > getRcp() 
    {return rcp(this);}
  //@}
    
protected:

  /** */
  virtual SolverState<Scalar> solveUnprec(const LinearOperator<Scalar>& op,
    const Vector<Scalar>& rhs,
    Vector<Scalar>& soln) const ;

    
};


}

#endif
