/* @HEADER@ */
//
 /* @HEADER@ */

#ifndef PLAYA_NONLINEAROPERATORBASE_HPP
#define PLAYA_NONLINEAROPERATORBASE_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaOut.hpp"
#include "PlayaObjectWithVerbosity.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearCombinationDecl.hpp"

namespace Playa
{
using namespace Teuchos;

/** 
 * Base class for nonlinear operators
 */
template <class Scalar>
class NonlinearOperatorBase 
  : public Handleable<NonlinearOperatorBase<Scalar> >,
    public ObjectWithVerbosity
{
public:
  /** Empty ctor, for contexts in which we don't know the
   * domain and range spaces at the beginning of construction time */
  NonlinearOperatorBase() 
    : domain_(), range_(), 
      jacobianIsValid_(false),
      residualIsValid_(false),
      currentEvalPt_(),
      currentFunctionValue_(),
      currentJ_()
    {;}

  /** Construct a nonlinear operator with a domain and range */
  NonlinearOperatorBase(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range) 
    : domain_(domain.ptr()), range_(range.ptr()), 
      jacobianIsValid_(false),
      residualIsValid_(false),
      currentEvalPt_(),
      currentFunctionValue_(),
      currentJ_()
    {;}
                            
  /** Return the domain space */
  const RCP<const VectorSpaceBase<Scalar> >& domain() const 
    {return domain_;}

  /** Return the range space */
  const RCP<const VectorSpaceBase<Scalar> >& range() const 
    {return range_;}

  /** Set the evaluation point */
  void setEvalPt(const Vector<Scalar>& x)
    {
      if (this->verb() >= 1)
      {
        Out::os() << "NonlinearOperatorBase Setting new eval pt";
        if (this->verb() > 3)
        {
          Out::os() << " to " << std::endl ;
          x.print(Out::os());
        }
        Out::os() << std::endl;
      }
      jacobianIsValid_ = false;
      residualIsValid_ = false;

      currentEvalPt_ = x.copy();
    }

  /** Get the current point at which the function is to be 
   * evaluated */
  const Vector<double>& currentEvalPt() const {return currentEvalPt_;}

  /** Return the Jacobian at the current evaluation point */
  LinearOperator<double> getJacobian() const 
    {
      if (this->verb() > 1)
      {
        Out::os() << "NonlinearOperatorBase getting Jacobian" << std::endl;
      }
      if (!jacobianIsValid_)
      {
        if (this->verb() > 3)
        {
          Out::os() << "...computing new J and F" << std::endl;
        }
        currentJ_ 
          = computeJacobianAndFunction(currentFunctionValue_);
        jacobianIsValid_ = true;
        residualIsValid_ = true;
      }
      else
      {
        if (this->verb() > 1)
        {
          Out::os() << "...reusing valid J" << std::endl;
        }
      }
      if (this->verb() > 3)
      {
        Out::os() << "J is " << std::endl;
        currentJ_.print(Out::os());
        Out::os() << std::endl;
      }
      return currentJ_;
    }

      

  /** Return the function value at the current evaluation point */
  Vector<double> getFunctionValue() const 
    {
      if (this->verb() > 1)
      {
        Out::os() << "NonlinearOperatorBase getting function value" << std::endl;
      }
      if (!residualIsValid_)
      {
        if (this->verb() > 1)
        {
          Out::os() << "...computing new F" << std::endl;
        }
        currentFunctionValue_ = computeFunctionValue();
        residualIsValid_ = true;
      }
      else
      {
        if (this->verb() > 1)
        {
          Out::os() << "...reusing valid F" << std::endl;
        }
      }

      if (this->verb() > 3)
      {
        Out::os() << "F is " << std::endl;
        currentFunctionValue_.print(Out::os());
        Out::os() << std::endl;
      }
      return currentFunctionValue_;
    }


  /** Return an initial guess appropriate to this problem */
  virtual Vector<double> getInitialGuess() const = 0 ;


protected:

  /** Compute the Jacobian at the current eval point */
  virtual LinearOperator<Scalar> computeJacobianAndFunction(Vector<double>& functionValue) const = 0 ;

  /** Compute the function value at the current eval point */
  virtual Vector<Scalar> computeFunctionValue() const 
    {
      computeJacobianAndFunction(currentFunctionValue_);
      return currentFunctionValue_;
    }

      
  /** Set the domain and range. This is protected so that solver
   * developers don't try to change the spaces on the fly */
  void setDomainAndRange(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range)
    {
      domain_ = domain.ptr();
      range_ = range.ptr();
    }


private:
  /** */
  RCP<const VectorSpaceBase<Scalar> > domain_;

  /** */
  RCP<const VectorSpaceBase<Scalar> > range_;

  /** */
  mutable bool jacobianIsValid_;

  /** */
  mutable bool residualIsValid_;

  /** */
  Vector<double> currentEvalPt_;

  /** */
  mutable Vector<double> currentFunctionValue_;

  /** */
  mutable LinearOperator<double> currentJ_;
};



 
}


#endif
