/* @HEADER@ */
//   
 /* @HEADER@ */


#ifndef POISSON_BOLTZMANN_JACOBIAN_HPP
#define POISSON_BOLTZMANN_JACOBIAN_HPP


#include "PlayaOperatorBuilder.hpp"

using namespace Playa;
using namespace Teuchos;


namespace Playa
{
  /** */
  class PoissonBoltzmannJacobian : public OperatorBuilder<double>
  {
  public:
    
    /** */
    PoissonBoltzmannJacobian(int nLocal, const VectorType<double>& vecType);

    /** */
    void setEvalPoint(const Vector<double>& x);

    /** */
    LinearOperator<double> getOp() const {return op_;}

    /** */
    const double& h() const {return h_;}

    /** */
    int nLocalRows() const {return nLocalRows_;}

  private:
    LinearOperator<double> op_;

    int nLocalRows_;

    double h_;
  };

}

#endif
