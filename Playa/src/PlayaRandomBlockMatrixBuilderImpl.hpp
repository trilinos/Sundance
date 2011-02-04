/* @HEADER@ */
//   
 /* @HEADER@ */


#ifndef RANDOMBLOCKMATRIX_BUILDER_IMPL_HPP
#define RANDOMBLOCKMATRIX_BUILDER_IMPL_HPP

#include "PlayaRandomBlockMatrixBuilderDecl.hpp"
#include "PlayaSimpleBlockOpDecl.hpp"
#include "PlayaRandomSparseMatrixBuilderDecl.hpp"



#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaSimpleBlockOpImpl.hpp"
#include "PlayaRandomSparseMatrixBuilderImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;


namespace Playa
{

template <class Scalar> 
inline RandomBlockMatrixBuilder<Scalar>
::RandomBlockMatrixBuilder(const VectorSpace<Scalar>& d,
  const VectorSpace<Scalar>& r,
  double blockDensity,
  double onProcDensity,
  double offProcDensity,
  const VectorType<double>& type)
  : OperatorBuilder<double>(d, r, type), op_()
{
  RCP<SimpleBlockOp<Scalar> > b = 
    rcp(new SimpleBlockOp<Scalar>(this->domain(), this->range()));
  RCP<LinearOperatorBase<Scalar> > p = b;
  op_ = p;

  for (int i=0; i<this->range().numBlocks(); i++)
  {
    for (int j=0; j<this->domain().numBlocks(); j++)
    {
      RandomSparseMatrixBuilder<Scalar> builder(this->domain().getBlock(j),
        this->range().getBlock(i),
        onProcDensity,
        offProcDensity,
        this->vecType());
      op_.setBlock(i, j, builder.getOp());
    }
  }
  b->endBlockFill();
}
}

#endif
