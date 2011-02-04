/* @HEADER@ */
//   
 /* @HEADER@ */


#ifndef RANDOMSPARSEMATRIX_BUILDER_DECL_HPP
#define RANDOMSPARSEMATRIX_BUILDER_DECL_HPP

#include "PlayaOperatorBuilder.hpp"


using namespace Playa;
using namespace Teuchos;


namespace Playa
{
/** */
template <class Scalar>
class RandomSparseMatrixBuilder : public OperatorBuilder<Scalar>
{
public:
  /** */
  RandomSparseMatrixBuilder(int nLocalRows, int nLocalCols,
    double onProcDensity,
    double offProcDensity,
    const VectorType<double>& vecType);
  /** */
  RandomSparseMatrixBuilder(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range,
    double onProcDensity,
    double offProcDensity,
    const VectorType<double>& vecType);

  /** */
  LinearOperator<double> getOp() const {return op_;}

private:

  void initOp(double onProcDensity,
    double offProcDensity);

  LinearOperator<double> op_;
};

}

#endif
