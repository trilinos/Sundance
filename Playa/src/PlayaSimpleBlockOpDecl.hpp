/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SIMPLEBLOCKOP_DECL_HPP
#define PLAYA_SIMPLEBLOCKOP_DECL_HPP


#include "PlayaDefs.hpp"
#include "Teuchos_Array.hpp"
#include "PlayaBlockOperatorBaseDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"



namespace Playa
{
using namespace Teuchos;

/**
 * Array-based block operator
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class SimpleBlockOp : public LinearOpWithSpaces<Scalar>,
                      public SetableBlockOperatorBase<Scalar>
{
public:
  /** */
  SimpleBlockOp(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range);

  /** */
  int numBlockRows() const  ;

  /** */
  int numBlockCols() const ;

  /** */
  const LinearOperator<Scalar>& getBlock(int i, int j) const ;

  /** */
  LinearOperator<Scalar> getNonconstBlock(int i, int j) ;

  /** */
  void setBlock(int i, int j, const LinearOperator<Scalar>& Aij) ; 

  /** */
  void apply(Teuchos::ETransp transApplyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const ;

private:
  Array<Array<LinearOperator<Scalar> > > blocks_;

}; 


/** \relates SimpleBlockOp Nonmember function to create a SimpleBlockOp */
template <class Scalar>
LinearOperator<Scalar> makeBlockOperator(
  const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range
  );





}

#endif
