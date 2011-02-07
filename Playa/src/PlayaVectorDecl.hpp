/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTORDECL_HPP
#define PLAYA_VECTORDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaVectorBaseDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaVectorFunctorsDecl.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace PlayaExprTemplates
{
template <class Scalar, class Node1, class Node2> class LC2;
template <class Scalar, class Node> class OpTimesLC; 
/** 
 * 
 */
enum LCSign {LCAdd = 1, LCSubtract = -1};
}



namespace Playa
{
  
/** 
 * User-level vector class. 
 *
 * <h2> Vector creation </h2>
 *
 * Ordinarily, you will never construct a Vector directly
 * from a derived type.  Rather, the createMember() method of
 * VectorSpace is used to build a vector of the appropriate
 * type, for example,
 * \code
 * VectorType<double> vecType = new EpetraVectorType();
 * int dimension = 100;
 * VectorSpace<double> space = vecType.createSpace(dimension);
 * Vector<double> x = space.createMember();
 * Vector<double> y = space.createMember();
 * \endcode
 * This hides from you all the ugly
 * details of creating a particular concrete type.
 *
 * You will frequently create an empty vector to be filled in later,
 * for example,
 * \code
 * Vector<double> y;
 * \endcode
 * Note that this vector isn't just empty, it's null. Not only does
 * it have no values assigned, it does not have a concrete type. An
 * call a method on a null vector will result in an error. What you
 * <it>can</it> do with a null vector is
 * <ul>
 * <li> assign another vector to it
 * \code
 * Vector<double> x = space.createVector();
 * Vector<Scalar> y;
 * y = x.copy();
 * \endcode
 * <li> assign the result of a vector operation to it
 * \code
 * Vector<Scalar> z = a*x + b*y;
 * \endcode
 *
 * <h2> Vector creation </h2>
 */
template <class Scalar>
class Vector : public Playa::Handle<VectorBase<Scalar> >
{
public:
  /** \name Constructors, Destructors, and Assignment Operators */
  //@{
  HANDLE_CTORS(Vector<Scalar>, VectorBase<Scalar>);

  /** Construct a vector from a 2-term LC */
  template<class Node1, class Node2>
  Vector(const PlayaExprTemplates::LC2<Scalar, Node1, Node2>& x);

  /** Construct a vector from an operator times a linear combination */
  template<class Node>
  Vector(const PlayaExprTemplates::OpTimesLC<Scalar, Node>& x);

  /** Assign a linear combination of vectors to this vector */
  template<class Node1, class Node2>
  Vector& operator=(const PlayaExprTemplates::LC2<Scalar, Node1, Node2>& x);

  /** Assign a scaled linear combination to this vector */
  template<class Node>
  Vector& operator=(const PlayaExprTemplates::OpTimesLC<Scalar, Node>& x);
  //@}

  /** \name Structural information */
  //@{
  /** My space */
  VectorSpace<Scalar> space() const 
    {return this->ptr()->space();}

  /** My communicator */
  const MPIComm& comm() const 
    {return this->space().comm();}

  /** My dimension  */
  int dim() const
    {
      return this->ptr()->space()->dim();
    }
  //@}
      

  /** \name Block operations */
  //@{
  /** get number of blocks */
  int numBlocks() const ;

  /** set the i-th block  */
  void setBlock(int i, const Vector<Scalar>& v);
      
  /** get the i-th block */
  const Vector<Scalar>& getBlock(int i) const;
      
  /** get the i-th block */
  Vector<Scalar> getNonConstBlock(int i) ;
      
  /** get the i-th block */
  const Vector<Scalar>& getBlock(const BlockIterator<Scalar>& b) const;
      
  /** get the i-th block */
  Vector<Scalar> getNonConstBlock(const BlockIterator<Scalar>& b);
  //@}



  /** \name Sequential data accessors */
  //@{
  /** Get the next chunk of values for read-only access */
  ConstDataChunk<Scalar> nextConstChunk() const ;
    
  /** Get the next chunk of values for read-write access */
  NonConstDataChunk<Scalar> nextChunk() ;

  /** Tell whether there are more chunks remaining to be accessed */
  bool hasMoreChunks() const ;

  /** Reset the data stream back to a state where all chunks are
   * considered unread. */
  void rewind() const ;
  //@}

  /** \name Random access to local elements */
  //@{
  /** const bracket operator for read-only random access to
   * local elements as specified by
   * a flat index that runs from 0 to space().numLocalElements(). 
   * If the vector does not consist of a single contiguous data chunk,
   * this might be slow (worst case would be O(N), if every element
   * is stored in its own chunk of length 1). 
   */
  const Scalar& operator[](int localIndex) const ;
    
  /** non-const bracket operator for read-write random access to
   * local elements as specified by
   * a flat index that runs from 0 to space().numLocalElements(). 
   * If the vector does not consist of a single contiguous data chunk,
   * this might be slow (worst case would be O(N), if every element
   * is stored in its own chunk of length 1). 
   */
  Scalar& operator[](int localIndex);

  /** parentheses operator for read-only random access to
   * local elements as specified by 
   * a block iterator and a flat index indicating the
   * element's location within that block. 
   */
  const Scalar& operator()(const BlockIterator<Scalar>& b,
    int localIndexWithinBlock) const ;
    
/** parentheses operator for read-write random access to
   * local elements as specified by 
   * a block iterator and a flat index indicating the
   * element's location within that block. 
   */
  Scalar& operator()(const BlockIterator<Scalar>& b,
    int localIndexWithinBlock) ;
  //@}

  /** \name Diagnostic output */
  //@{
  /** Return a short string description of the vector */
  std::string description() const ;

  /** Print the vector in some detail */
  void print(std::ostream& os) const ;    
  //@}


  /** \name Math operations */
  //@{
  /** Multiply this vector by a constant scalar factor 
   * \code
   * this = alpha * this;
   * \endcode
   */
  Vector<Scalar>& scale(const Scalar& alpha);

  /** 
   * Add a scaled vector to this vector times a constant:
   * \code
   * this = gamma*this + alpha*x 
   * \endcode
   */
  Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
    const Scalar& gamma=1.0);
  /** 
   * Add two scaled vectors to this vector times a constant:
   * \code
   * this = alpha*x + beta*y + gamma*this
   * \endcode
   */
  Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
    const Scalar& beta, const Vector<Scalar>& y, 
    const Scalar& gamma);

  /** 
   * Copy the values of another vector into this vector
   * \code
   * this = x
   * \endcode
   */
  Vector<Scalar>& acceptCopyOf(const Vector<Scalar>& x);

  /** 
   * Create a new vector that is a copy of this vector 
   */
  Vector<Scalar> copy() const ;

  /** 
   * In-place element-by-element product (Matlab dot-star operator)
   */
  Vector<Scalar>& dotStar(const Vector<Scalar>& other) ;

  /** 
   * In-place element-by-element division (Matlab dot-slash operator)
   */
  Vector<Scalar>& dotSlash(const Vector<Scalar>& other) ;

  /** 
   * Element-by-element product (Matlab dot-star operator)
   */
  Vector<Scalar> dotStar(const Vector<Scalar>& other) const ;

  /** 
   * Element-by-element division (Matlab dot-slash operator)
   */
  Vector<Scalar> dotSlash(const Vector<Scalar>& other) const ;

  /** 
   * Return element-by-element reciprocal as a new vector
   */
  Vector<Scalar> reciprocal() const ;


  /** 
   * Return element-by-element absolute value as a new vector
   */
  Vector<Scalar> abs() const ;

  /** 
   * Overwrite self with element-by-element reciprocal
   */
  Vector<Scalar>& reciprocal() ;

  /** 
   * Overwrite self with element-by-element absolute value 
   */
  Vector<Scalar>& abs() ;

  /** 
   * Overwrite self with random values
   */
  Vector<Scalar>& randomize() ;

  /** 
   * Set all elements to a constant value
   */
  void setToConstant(const Scalar& alpha) ;

      
  /** 
   * Take dot product with another vector
   */
  Scalar dot(const Vector<Scalar>& other) const ;

  /** 
   * Overloaded operator for dot product 
   */
  Scalar operator*(const Vector<Scalar>& other) const ;

  /**
   * Compute the 1-norm of this vector
   */
  Scalar norm1() const ;

  /**
   * Compute the 2-norm of this vector
   */
  Scalar norm2() const ;

  /**
   * Compute the weighted 2-norm of this vector
   */
  Scalar norm2(const Vector<Scalar>& weights) const ;    


  /**
   * Compute the infinity-norm of this vector
   */
  Scalar normInf() const ;

  /**
   * Set all elements to zero 
   */
  void zero();


  /** Return the max element */
  Scalar max() const;

  /** Return the min element */
  Scalar min()const;


  //@}




  /** Get a stopwtach for timing vector operations */
  static RCP<Time>& opTimer()
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("Low-level vector operations");
      return rtn;
    }

  Vector<Scalar> eval() const {return copy();}

  bool containsVector(const VectorBase<Scalar>* vec) const
    {return this->ptr().get()==vec;}

  void evalInto(Vector<Scalar>& other) const {other.acceptCopyOf(*this);}

  void addInto(Vector<Scalar>& other, PlayaExprTemplates::LCSign sign) const
    {
      other.update(sign, *this);
    }

  /** \name Functor application */
  //@{

  /** Apply a unary functor, overwriting this vector with the results */
  template <class UF> 
  Vector<Scalar>& applyUnaryFunctor(const UF& functor);

  /** Apply a unary functor to another vector, writing the results
      into this vector. The other vector is unchanged. */
  template <class UF> 
  Vector<Scalar>& acceptUnaryFunctor(const UF& functor,
    const Vector<Scalar>& other);

  /** Apply a binary functor to this vector and another vector, 
      writing the results
      into this vector. The other vector is unchanged. */
  template <class VF> 
  Vector<Scalar>& applyBinaryFunctor(const VF& functor,
    const Vector<Scalar>& other);

  /** Apply a ternary functor to this vector and two other vectors, 
      writing the results
      into this vector. The other vectors are unchanged. */
  template <class VF> 
  Vector<Scalar>& applyTernaryFunctor(const VF& functor,
    const Vector<Scalar>& x, const Vector<Scalar>& y);

  /** Apply a unary reduction functor */
  template <class RF> 
  typename PlayaFunctors::VectorFunctorTraits<Scalar, RF>::ReturnType 
  applyUnaryReductionFunctor(
    const RF& functor)
    const ;

  /** Apply a binary reduction functor */
  template <class RF> 
  typename PlayaFunctors::VectorFunctorTraits<Scalar, RF>::ReturnType 
  applyBinaryReductionFunctor(const RF& functor, const Vector<Scalar>& other)
    const ;
    

    
  //@}
    
protected:

  /** get a subblock as specified by a deque of indices */
  const Vector<Scalar>& getBlock(const std::deque<int>& iter) const ;

  /** get a non-const subblock as specified by a deque of indices */
  Vector<Scalar> getNonConstBlock(const std::deque<int>& iter) ;
  

private:

};


template <class Scalar> class LoadableVector;
/** \relates Vector */
template <class Scalar>
LoadableVector<Scalar>* loadable(Vector<Scalar> vec) ;


}

template <class Scalar> inline
std::ostream& operator<<(std::ostream& os, const Playa::Vector<Scalar>& x) 
{
  x.print(os);
  return os;
}




#endif
