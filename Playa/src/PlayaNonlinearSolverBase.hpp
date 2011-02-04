/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_NONLINEARSOLVERBASE_HPP
#define PLAYA_NONLINEARSOLVERBASE_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"

namespace Playa
{
  using namespace Teuchos;

  /*! \brief Here I am
   *
   */
  template <class Scalar>
  class NonlinearSolverBase 
  {
  public:
    /** */
    NonlinearSolverBase(const ParameterList& params = ParameterList());

    /** */
    virtual ~NonlinearSolverBase(){;}

    /** */
    virtual Vector<double> solve() const = 0  ;

  protected:

    const ParameterList& params() const {return params_;}

  private:
    ParameterList params_;
  };

  
  template <class Scalar> inline
  NonlinearSolverBase<Scalar>
  ::NonlinearSolverBase(const ParameterList& params)
    : params_(params)
  {;}
  
}

#endif
