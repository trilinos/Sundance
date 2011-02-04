/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_LINEARCOMBINATIONDECL_HPP
#define PLAYA_LINEARCOMBINATIONDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"



namespace PlayaOps
{
using Playa::Vector;
using Playa::VectorBase;
using Playa::VectorSpace;
using Playa::LinearOperator;
/** 
 *
 */
template <class Scalar> 
class ConvertibleToVector
{
public:
  /** */
  virtual ~ConvertibleToVector(){;}

  /** */
  virtual Vector<Scalar> eval() const = 0 ;

  /** */
  VectorSpace<Scalar> space() const {return eval().space();}

  /** Return the dimension of the vector  */
  int dim() const {return eval().dim();} 

  /** 
   * Create a new vector that is a copy of this vector 
   */
  Vector<Scalar> copy() const {return eval().copy();}

  /** 
   * Element-by-element product (Matlab dot-star operator)
   */
  Vector<Scalar> dotStar(const Vector<Scalar>& other) const 
    {return eval().dotStar(other);}

  /** 
   * Element-by-element division (Matlab dot-slash operator)
   */
  Vector<Scalar> dotSlash(const Vector<Scalar>& other) const 
    {return eval().dotSlash(other);}

  /** 
   * Return element-by-element reciprocal as a new vector
   */
  Vector<Scalar> reciprocal() const {return reciprocal();}

  /** 
   * Return element-by-element absolute value as a new vector
   */
  Vector<Scalar> abs() const {return abs();} 

  /** */
  Scalar norm2() const {return eval().norm2();}

  /** */
  Scalar norm1() const {return eval().norm1();}

  /** */
  Scalar normInf() const {return eval().normInf();}

  /** */
  Scalar max() const {return eval().max();}

  /** */
  Scalar min() const {return eval().min();}

  /** Return the min element and the corresponding index */
  Scalar min(int& index)const {return eval().min(index);}

  /** Return the max element and the corresponding index */
  Scalar max(int& index)const {return eval().max(index);}
};

/** 
 * Class OpTimesLC holds an operator times something convertible to a vector
 */
template <class Scalar, class Node>
class OpTimesLC : public ConvertibleToVector<Scalar>
{
public:

  /** */
  virtual ~OpTimesLC(){;}

  /** */
  OpTimesLC(const Scalar& alpha, const Node& x);

  /** */
  OpTimesLC(const Scalar& alpha,
    const Playa::LinearOperator<Scalar>& op, 
    const Node& x);

  /** 
   * Evaluate the term into the argument vector, overwriting 
   * the previous value of the argument. */
  void evalInto(Playa::Vector<Scalar>& result) const ;

  /** Add the term into the argument vector */
  void addInto(Playa::Vector<Scalar>& result, 
    LCSign sign = LCAdd) const ;

  /** Evaluate the term and return its value */
  virtual Playa::Vector<Scalar> eval() const ;

  /** Determine whether this term contains the given vector */
  bool containsVector(const VectorBase<Scalar>* vec) const ;

  /** */
  const LinearOperator<Scalar>& op() const {return op_;}

  /** */
  const Scalar& alpha() const {return alpha_;}

  /** */
  const Node& node() const {return x_;}

  /** */
  Scalar norm2() const {return eval().norm2();}

  /** */
  Scalar normInf() const {return eval().normInf();}
    
private:
  Scalar alpha_;
    
  LinearOperator<Scalar> op_;

  Node x_;

  /** */
  static Scalar one() {return Teuchos::ScalarTraits<Scalar>::one();}

  /** */
  static Scalar zero() {return Teuchos::ScalarTraits<Scalar>::zero();}
};


/**
 * Class LC2 is a 2-term linear combination
 */
template <class Scalar, class Node1, class Node2>
class LC2  : public ConvertibleToVector<Scalar>
{
public:
  /** */
  virtual ~LC2(){;}

  /** */
  LC2(const Node1& x1, const Node2& x2, LCSign sign = LCAdd);

  /** */
  void evalInto(Playa::Vector<Scalar>& result) const ;

  /** */
  void addInto(Playa::Vector<Scalar>& result, 
    LCSign sign = LCAdd) const ;

  /** */
  virtual Playa::Vector<Scalar> eval() const ;

  /** */
  bool containsVector(const VectorBase<Scalar>* vec) const ;

    
private:
  Node1 x1_;

  Node2 x2_;

  LCSign sign_;

  /** */
  static Scalar one() {return Teuchos::ScalarTraits<Scalar>::one();}

  /** */
  static Scalar zero() {return Teuchos::ScalarTraits<Scalar>::zero();}
};


}

namespace Playa
{
using PlayaOps::OpTimesLC;
using PlayaOps::LC2;
using PlayaOps::LCAdd;
using PlayaOps::LCSubtract;

/* ------------------------ global methods ----------------------- */


/*======================================================================
 *
 *    scalar times vector
 *
 *======================================================================*/

/* scalar * vec */
template <class Scalar> 
OpTimesLC<Scalar, Vector<Scalar> > operator*(const Scalar& alpha, 
  const Vector<Scalar>& x);

/* vec * scalar */
template <class Scalar> 
OpTimesLC<Scalar, Vector<Scalar> > operator*(const Vector<Scalar>& x, 
  const Scalar& alpha);


/*======================================================================
 *
 *    scalar times OpTimesLC
 *
 *======================================================================*/

/* scalar * OpTimesLC */
template <class Scalar, class Node> 
OpTimesLC<Scalar, Node> 
operator*(const Scalar& alpha, 
  const OpTimesLC<Scalar, Node>& x);

/* OpTimesLC * scalar */
template <class Scalar, class Node> 
OpTimesLC<Scalar, Node> 
operator*(const OpTimesLC<Scalar, Node>& x, const Scalar& alpha);


/*======================================================================
 *
 *    scalar times LC2
 *
 *======================================================================*/

/* scalar * LC2 */
template <class Scalar, class Node1, class Node2> 
OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> > 
operator*(const Scalar& alpha, 
  const LC2<Scalar, Node1, Node2>& x);

/* LC2 * scalar */
template <class Scalar, class Node1, class Node2> 
OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> > 
operator*(const LC2<Scalar, Node1, Node2>& x, const Scalar& alpha);
  


/*======================================================================
 *
 *    operator times [vectors, OpTimesLC, LC2]
 *
 *======================================================================*/

/* op * vec */
template <class Scalar>
OpTimesLC<Scalar, Vector<Scalar> > 
operator*(const LinearOperator<Scalar>& op, 
  const Vector<Scalar>& x);


/* op * OpTimesLC */
template <class Scalar, class Node> 
OpTimesLC<Scalar, Node> 
operator*(const LinearOperator<Scalar>& op, 
  const OpTimesLC<Scalar, Node>& x);


/* op * LC2 */
template <class Scalar, class Node1, class Node2> 
OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> > 
operator*(const LinearOperator<Scalar>& op, 
  const LC2<Scalar, Node1, Node2>& x);


/*======================================================================
 *
 *    add/subtract vector, vector
 *
 *======================================================================*/
  
/* vec + vec */
template <class Scalar> 
LC2<Scalar, Vector<Scalar>, Vector<Scalar> >
operator+(const Vector<Scalar>& x1, 
  const Vector<Scalar>& x2);
  
/* vec - vec */
template <class Scalar> 
LC2<Scalar, Vector<Scalar>, Vector<Scalar> >
operator-(const Vector<Scalar>& x1, 
  const Vector<Scalar>& x2);


/*======================================================================
 *
 *    add/subtract vector, OpTimesLC
 *
 *======================================================================*/


/* vec + OpTimesLC */
template <class Scalar, class Node> 
LC2<Scalar, Vector<Scalar>, OpTimesLC<Scalar, Node> >
operator+(const Vector<Scalar>& x1, 
  const OpTimesLC<Scalar, Node>& x2);


/* vec - OpTimesLC */
template <class Scalar, class Node> 
LC2<Scalar, Vector<Scalar>, OpTimesLC<Scalar, Node> >
operator-(const Vector<Scalar>& x1, 
  const OpTimesLC<Scalar, Node>& x2);


/* OpTimesLC + vec */
template <class Scalar, class Node> 
LC2<Scalar, OpTimesLC<Scalar, Node>, Vector<Scalar> >
operator+(const OpTimesLC<Scalar, Node>& x1, 
  const Vector<Scalar>& x2);
  
/* OpTimesLC - vec */
template <class Scalar, class Node> 
LC2<Scalar, OpTimesLC<Scalar, Node>, Vector<Scalar> >
operator-(const OpTimesLC<Scalar, Node>& x1, 
  const Vector<Scalar>& x2);


  
/*======================================================================
 *
 *    add/subtract OpTimesLC, OpTimesLC
 *
 *======================================================================*/
  
/* OpTimesLC + OpTimesLC */
template <class Scalar, class Node1, class Node2> 
LC2<Scalar, OpTimesLC<Scalar, Node1>, OpTimesLC<Scalar, Node2> >
operator+(const OpTimesLC<Scalar, Node1>& x1, 
  const OpTimesLC<Scalar, Node2>& x2);

  
/* OpTimesLC - OpTimesLC */
template <class Scalar, class Node1, class Node2> 
LC2<Scalar, OpTimesLC<Scalar, Node1>, OpTimesLC<Scalar, Node2> >
operator-(const OpTimesLC<Scalar, Node1>& x1, 
  const OpTimesLC<Scalar, Node2>& x2);
  

  
/*======================================================================
 *
 *    add/subtract Vector, LC2
 *
 *======================================================================*/

  
/* vec + LC2 */
template <class Scalar, class Node1, class Node2> 
LC2<Scalar, Vector<Scalar>, LC2<Scalar, Node1, Node2> >
operator+(const Vector<Scalar>& x1, 
  const LC2<Scalar, Node1, Node2>& x2);


/* vec - LC2 */
template <class Scalar, class Node1, class Node2> 
LC2<Scalar, Vector<Scalar>, LC2<Scalar, Node1, Node2> >
operator-(const Vector<Scalar>& x1, 
  const LC2<Scalar, Node1, Node2>& x2);


/* LC2 + vec */
template <class Scalar, class Node1, class Node2> 
LC2<Scalar, LC2<Scalar, Node1, Node2>, Vector<Scalar> >
operator+(const LC2<Scalar, Node1, Node2>& x1, 
  const Vector<Scalar>& x2);


/* LC2 - vec */
template <class Scalar, class Node1, class Node2> 
LC2<Scalar, LC2<Scalar, Node1, Node2>, Vector<Scalar> >
operator-(const LC2<Scalar, Node1, Node2>& x1, 
  const Vector<Scalar>& x2);



/*======================================================================
 *
 *    add/subtract OpTimesLC, LC2
 *
 *======================================================================*/


/* OpTimesLC + LC2 */
template <class Scalar, class Node0, class Node1, class Node2> 
LC2<Scalar, OpTimesLC<Scalar, Node0>, LC2<Scalar, Node1, Node2> > 
operator+(const OpTimesLC<Scalar, Node0>& x1, 
  const LC2<Scalar, Node1, Node2>& x2);

/* OpTimesLC - LC2 */
template <class Scalar, class Node0, class Node1, class Node2> 
LC2<Scalar, OpTimesLC<Scalar, Node0>, LC2<Scalar, Node1, Node2> > 
operator-(const OpTimesLC<Scalar, Node0>& x1, 
  const LC2<Scalar, Node1, Node2>& x2);



/* LC2 + OpTimesLC */
template <class Scalar, class Node1, class Node2, class Node3> 
LC2<Scalar, LC2<Scalar, Node1, Node2>, OpTimesLC<Scalar, Node3> > 
operator+(const LC2<Scalar, Node1, Node2>& x1, 
  const OpTimesLC<Scalar, Node3>& x2);


/* LC2 - OpTimesLC */
template <class Scalar, class Node1, class Node2, class Node3> 
LC2<Scalar, LC2<Scalar, Node1, Node2>, OpTimesLC<Scalar, Node3> > 
operator-(const LC2<Scalar, Node1, Node2>& x1, 
  const OpTimesLC<Scalar, Node3>& x2);



/*======================================================================
 *
 *    add/subtract LC2, LC2
 *
 *======================================================================*/
  
/* LC2 + LC2 */
template <class Scalar, class Node1, class Node2, 
          class Node3, class Node4> 
LC2<Scalar, LC2<Scalar, Node1, Node2>, LC2<Scalar, Node3, Node4> >
operator+(const LC2<Scalar, Node1, Node2>& x1, 
  const LC2<Scalar, Node3, Node4>& x2);


/* LC2 - LC2 */
template <class Scalar, class Node1, class Node2, 
          class Node3, class Node4> 
LC2<Scalar, LC2<Scalar, Node1, Node2>, LC2<Scalar, Node3, Node4> >
operator-(const LC2<Scalar, Node1, Node2>& x1, 
  const LC2<Scalar, Node3, Node4>& x2);



  


  


}



#endif
