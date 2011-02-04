/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_LINEAROPERATORDECL_HPP
#define PLAYA_LINEAROPERATORDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaLinearOperatorBaseDecl.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "PlayaRowAccessibleOp.hpp"


namespace Playa
{
using namespace Teuchos;


template <class Scalar>  class LinearSolver;
template <class Scalar>  class VectorSpace;
template <class Scalar>  class Vector;
template <class Scalar>  class VectorType;

/** 
 * User-level linear operator class
 */
template <class Scalar>
class LinearOperator : public Playa::Handle<LinearOperatorBase<Scalar> >
{
public:
  /** \name Constructors, Destructors, and Assignment Operators */
  //@{
  /** Empty constructor*/
  LinearOperator();

  /** Constructor with smart pointer */
  LinearOperator(const RCP<LinearOperatorBase<Scalar> >& smartPtr);
  //@}

  /** Return the domain */
  const VectorSpace<Scalar> domain() const ;

  /** Return the range */
  const VectorSpace<Scalar> range() const ;


  /** 
   * Compute
   * \code
   * out = beta*out + alpha*op*in;
   * \endcode
   **/
  void apply(const Vector<Scalar>& in,
    Vector<Scalar>& out) const ;

  /**  
   * Compute
   * \code
   * out = beta*out + alpha*op^T*in;
   * \endcode
   **/
  void applyTranspose(const Vector<Scalar>& in,
    Vector<Scalar>& out) const ;


  //       /** For the moment this does nothing*/
  LinearOperator<Scalar> form() const {return *this;}
      
      
  /** Get a stopwatch for timing vector operations */
  RCP<Time>& opTimer();

  /**
   * Return a TransposeOperator.
   */
  LinearOperator<Scalar> transpose() const ; 


  /** Return a Loadable Matrix  */
  RCP<LoadableMatrix<Scalar> > matrix();

  /** Get a row of the underlying matrix */     
  void getRow(const int& row, 
    Teuchos::Array<int>& indices, 
    Teuchos::Array<Scalar>& values) const ;
    

  /** \name  Block operations  */
  //@{
      
  /** return number of block rows */
  int numBlockRows() const ;
      

  /** return number of block cols */
  int numBlockCols() const ;
      

  /** get the (i,j)-th block */
  LinearOperator<Scalar> getBlock(const int &i, const int &j) const ;


  /** get a writeable copy of the (i,j)-th block */
  LinearOperator<Scalar> getNonconstBlock(const int &i, const int &j) ;

  /** set the (i,j)-th block 
   *  If the domain and/or the range are not set, then we
   *  are building the operator
   */
  void setBlock(int i, int j, 
    const LinearOperator<Scalar>& sub);

  /** Finalize block assembly */
  void endBlockFill();

  //@}

      

private:

};

}


#endif
