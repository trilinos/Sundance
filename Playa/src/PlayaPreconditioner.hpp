/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_PRECONDITIONER_HPP
#define PLAYA_PRECONDITIONER_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaPreconditionerBase.hpp"

namespace Playa
{
  /**
   * 
   */
  template <class Scalar> 
  class Preconditioner : public Playa::Handle<PreconditionerBase<Scalar> >
  {
  public:
    /* Boilerplate ctors */
    HANDLE_CTORS(Preconditioner, PreconditionerBase<Scalar>);

    /** Change the value of a double parameter */
    void changeParameter(const std::string& name, const double& value);

    /** Change the value of an integer parameter */
    void changeParameter(const std::string& name, int value);

    
    
    /** Left preconditioner */
    LinearOperator<Scalar> left() const ;
    
    /** Right preconditioner */
    LinearOperator<Scalar> right() const ;
    
    /** return true if this preconditioner has both left and
     * right components. */
    bool isTwoSided() const {return hasLeft() && hasRight();}
    
    /** return true if this preconditioner has a nontrivial left component */
    bool hasLeft() const ;
    
    /** return true if this preconditioner has
     * a nontrivial right component */
    bool hasRight() const ;
    
    /** return true if this preconditioner has neither left nor
     * right operators defined */
    bool isIdentity() const {return !hasLeft() && !hasRight();}
  };

  

  template <class Scalar> inline 
  LinearOperator<Scalar> Preconditioner<Scalar>::left() const 
  {
    TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
                       "null pointer in Preconditioner<Scalar>::left()");
    return this->ptr()->left();
  }

  template <class Scalar> inline 
  LinearOperator<Scalar> Preconditioner<Scalar>::right() const 
  {
    TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
                       "null pointer in Preconditioner<Scalar>::right()");
    return this->ptr()->right();
  }

  template <class Scalar> inline
  bool Preconditioner<Scalar>::hasLeft() const 
  {
    return (this->ptr().get()!=0 && this->ptr()->hasLeft());
  }

  template <class Scalar> inline
  bool Preconditioner<Scalar>::hasRight() const 
  {
    return (this->ptr().get()!=0 && this->ptr()->hasRight());
  }

  
}

#endif
