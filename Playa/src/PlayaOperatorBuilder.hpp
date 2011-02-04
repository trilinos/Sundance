/* @HEADER@ */
//   
 /* @HEADER@ */


#ifndef PLAYA_OPERATORBUILDER_HPP
#define PLAYA_OPERATORBUILDER_HPP

#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearCombinationDecl.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_GlobalMPISession.hpp"


using namespace Playa;
using namespace Teuchos;


namespace Playa
{
  /** Base class for building test operators */
  template <class Scalar>
  class OperatorBuilder
  {
  public:
    /** */
    OperatorBuilder(int nLocal, const VectorType<Scalar>& vecType);
    /** */
    OperatorBuilder(int nLocalDomain, int nLocalRange,
                    const VectorType<Scalar>& vecType);
    /** */
    OperatorBuilder(const VectorSpace<Scalar>& domain,
                    const VectorSpace<Scalar>& range,
                    const VectorType<Scalar>& vecType);
    /** */
    virtual ~OperatorBuilder(){;}

    /** */
    const VectorType<Scalar>& vecType() const {return vecType_;}

    /** */
    const VectorSpace<Scalar>& domain() const {return domain_;}

    /** */
    const VectorSpace<Scalar>& range() const {return range_;}

    /** */
    virtual LinearOperator<Scalar> getOp() const = 0 ; 

  protected:

  private:
    VectorType<Scalar> vecType_;

    VectorSpace<Scalar> domain_;

    VectorSpace<Scalar> range_;
  };

  template <class Scalar> 
  inline OperatorBuilder<Scalar>
  ::OperatorBuilder(int nLocalRows, const VectorType<Scalar>& vecType)
    : vecType_(vecType), domain_(), range_()
  {
    range_ = vecType_.createEvenlyPartitionedSpace(MPIComm::world(), nLocalRows);
    domain_ = range_;
  }

  template <class Scalar> 
  inline OperatorBuilder<Scalar>
  ::OperatorBuilder(int nLocalDomain, int nLocalRange,
                    const VectorType<Scalar>& vecType)
    : vecType_(vecType), domain_(), range_()
  {
    range_ = vecType_.createEvenlyPartitionedSpace(MPIComm::world(), nLocalRange);
    domain_ = vecType_.createEvenlyPartitionedSpace(MPIComm::world(), nLocalDomain);
  }

  

  template <class Scalar> 
  inline OperatorBuilder<Scalar>
  ::OperatorBuilder(const VectorSpace<Scalar>& domain,
                    const VectorSpace<Scalar>& range,
                    const VectorType<Scalar>& vecType)
    : vecType_(vecType), domain_(domain), range_(range)
  {}

}

#endif
