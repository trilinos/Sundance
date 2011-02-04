/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_POISSONBOLTZMANNOP_HPP
#define PLAYA_POISSONBOLTZMANNOP_HPP

#include "PlayaDefs.hpp"
#include "PlayaNonlinearOperatorBase.hpp"
#include "PlayaPoissonBoltzmannJacobian.hpp"
#include "PlayaVectorDecl.hpp"

namespace Playa
{
  /** */
  class PoissonBoltzmannOp : public NonlinearOperatorBase<double>
  {
  public:
    /** */
    PoissonBoltzmannOp(int nLocal, const VectorType<double>& vecType);

    /** */
    Vector<double> getInitialGuess() const ;

    /** */
    Vector<double> exactSoln() const ;

    

    /* */
    GET_RCP(NonlinearOperatorBase<double>);

  protected:
    /** */
    LinearOperator<double> 
    computeJacobianAndFunction(Vector<double>& functionValue) const ;

  private:

    mutable PoissonBoltzmannJacobian J_;
    
    RCP<GhostImporter<double> > importer_;

    double uLeftBC_;

    double uRightBC_;

    
  };
 
}


#endif
