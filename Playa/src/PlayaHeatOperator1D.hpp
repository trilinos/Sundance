/* @HEADER@ */
//   
 /* @HEADER@ */


#ifndef PLAYA_HEATOPERATOR1D_HPP
#define PLAYA_HEATOPERATOR1D_HPP

#include "PlayaOperatorBuilder.hpp"

using namespace Playa;
using namespace Teuchos;


namespace Playa
{
  /** */
  class HeatOperator1D : public OperatorBuilder<double>
  {
  public:
    /** */
    HeatOperator1D(int nLocal, const VectorType<double>& vecType);

    /** */
    void set_cj(const double& cj) {cj_ = cj;}

    /** */
    LinearOperator<double> getOp() const ;

  private:
    RCP<MatrixFactory<double> > matrixFactory_;
    double cj_;
    double h_;
    int nLocalRows_;
    int lowestLocalRow_;
  };

}

#endif
