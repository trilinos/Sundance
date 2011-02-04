/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_ANASAZIEIGENSOLVER_DECL_HPP
#define PLAYA_ANASAZIEIGENSOLVER_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp" 
#include "PlayaSolverState.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaEigensolverBase.hpp"

namespace Playa
{
using Teuchos::ParameterList;

/**
 * Object wrapper for Anasazi eigenvalue solver.
 */
template <class Scalar>
class AnasaziEigensolver
  : public EigensolverBase<Scalar>,
    public Playa::Handleable<EigensolverBase<Scalar> >
{
public:
  /** */
  AnasaziEigensolver(const ParameterList& params) 
    : EigensolverBase<Scalar>(params) {;}

  /**
   * Solve a generalized eigenvalue problem \f$ K x = \lambda M x \f$
   */
  virtual void solve(
    const LinearOperator<Scalar>& K,
    const LinearOperator<Scalar>& M,
    Array<Vector<Scalar> >& ev,
    Array<std::complex<Scalar> >& ew) const ;

  /** \name Handleable interface */
  //@{
  /** Return a ref counted pointer to a newly created object */
  virtual RCP<EigensolverBase<Scalar> > getRcp() 
    {return rcp(this);}
  //@}
  
private:

  static Time& solveTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("AnasaziEigensolver solve()"); 
      return *rtn;
    }

  static Time& precondBuildTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("AnasaziEigensolver building preconditioner"); 
      return *rtn;
    }

};




}


#endif
