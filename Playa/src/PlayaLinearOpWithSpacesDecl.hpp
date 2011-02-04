/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_LINEAROPWITHSPACES_DECL_HPP
#define PLAYA_LINEAROPWITHSPACES_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearOperatorBaseDecl.hpp"

namespace Playa
{

/** 
 * LinearOpWithSpaces provides a simple implementation of the domain()
 * and range() methods of LinearOperatorBase.
 */
template <class Scalar>
class LinearOpWithSpaces : public LinearOperatorBase<Scalar> 
{
public:
  /** */
  LinearOpWithSpaces(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range)
    : domain_(domain.ptr()), range_(range.ptr()) {}

  /** Virtual dtor */
  ~LinearOpWithSpaces(){}

  /** Return the domain */
  const RCP<const VectorSpaceBase<Scalar> > domain() const 
    {return domain_;}

  /** Return the range */
  const RCP<const VectorSpaceBase<Scalar> > range() const 
    {return range_;}

private:
  RCP<const VectorSpaceBase<Scalar> > domain_;
  RCP<const VectorSpaceBase<Scalar> > range_;
};




}


#endif
