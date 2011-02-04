/* @HEADER@ */
//   
 /* @HEADER@ */


#ifndef RANDOMBLOCKMATRIX_BUILDER_DECL_HPP
#define RANDOMBLOCKMATRIX_BUILDER_DECL_HPP

#include "PlayaOperatorBuilder.hpp"
#include "PlayaRandomSparseMatrixBuilderDecl.hpp"

using namespace Playa;
using namespace Teuchos;


namespace Playa
{
  /** */
  template <class Scalar>
  class RandomBlockMatrixBuilder : public OperatorBuilder<Scalar>
  {
  public:
    /** */
    RandomBlockMatrixBuilder(const VectorSpace<Scalar>& domain,
      const VectorSpace<Scalar>& range,
      double blockDensity,
      double onProcDensity,
      double offProcDensity,
      const VectorType<double>& vecType);

    /** */
    LinearOperator<double> getOp() const {return op_;}

  private:
    LinearOperator<double> op_;
  };
}

#endif
