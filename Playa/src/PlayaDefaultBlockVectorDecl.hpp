/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_DEFAULT_BLOCK_VECTOR_DECL_HPP
#define PLAYA_DEFAULT_BLOCK_VECTOR_DECL_HPP

#include "PlayaBlockVectorBaseDecl.hpp"
#include "PlayaVectorDecl.hpp"
#include "Teuchos_Array.hpp"


namespace Playa
{
using Teuchos::Array;

/** 
 * Base class for blocked vectors 
 */
template <class Scalar>
class DefaultBlockVector : public BlockVectorBase<Scalar>
{
public:
  /** */
  DefaultBlockVector(const VectorSpace<Scalar>& space);

  /** */
  virtual ~DefaultBlockVector(){}

  /** \name VectorBase interface */
  //@{
  /** Access to the space in which this vector lives */
  RCP<const VectorSpaceBase<double> > space() const {return space_.ptr();}
  //@}

  /** */
  virtual void setBlock(int b, const Vector<Scalar>& block) ;

  /** */
  virtual const Vector<Scalar>& getBlock(int b) const ;

  /** */
  virtual Vector<Scalar> getNonConstBlock(int b) ;

  /** */
  virtual int numBlocks() const {return blocks_.size();}
  
private:
  VectorSpace<Scalar> space_;
  Teuchos::Array<Vector<Scalar> > blocks_;
  
};

}

#endif
