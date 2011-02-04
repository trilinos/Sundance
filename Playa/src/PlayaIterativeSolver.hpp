/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_ITERATIVESOLVER_HPP
#define PLAYA_ITERATIVESOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"

namespace Playa
{
  using namespace Teuchos;

  /**
   *
   */
  template <class Scalar>
  class IterativeSolver : public LinearSolverBase<Scalar>
  {
  public:
    /** */
    IterativeSolver(const ParameterList& params = ParameterList());

    /** */
    virtual ~IterativeSolver(){;}
    
    /** */
    int getMaxiters() const 
    {return this->parameters().template get<int>(maxitersParam());}

    /** */
    Scalar getTol() const 
    {return this->parameters().template get<double>(tolParam());}

    /** Change the convergence tolerance. */
    virtual void updateTolerance(const double& tol)
    {getParameter<double>(this->parameters(), tolParam())=tol;}

    /** */
    static std::string maxitersParam() {return "Max Iterations";}

    /** */
    static std::string tolParam() {return "Tolerance";}

    /** */
    static int defaultMaxiters() {return 500;}

    /** */
    static Scalar defaultTol() {return 1.0e-10;}
  };

  
  template <class Scalar> inline
  IterativeSolver<Scalar>::IterativeSolver(const ParameterList& params)
    : LinearSolverBase<Scalar>(params)
  {;}
  
}

#endif
