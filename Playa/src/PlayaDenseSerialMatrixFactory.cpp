/* @HEADER@ */
//   
 /* @HEADER@ */

#include "PlayaDenseSerialMatrixFactory.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaSerialVector.hpp"
#include "PlayaVectorSpaceDecl.hpp"  
#include "PlayaVectorDecl.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_MPIComm.hpp"
#include "PlayaLinearOperatorDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaVectorImpl.hpp"
#endif


namespace Playa
{

using namespace Teuchos;


DenseSerialMatrixFactory::DenseSerialMatrixFactory(
  const VectorSpace<double>& domain,
  const VectorSpace<double>& range)
  : 
  domain_(domain),
  range_(range)
{}


LinearOperator<double> DenseSerialMatrixFactory::createMatrix() const
{
  RCP<LinearOperatorBase<double> > A 
    = rcp(new DenseSerialMatrix(domain_, range_));
  return A;
}

}
