/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_ILUKPRECONDITIONERFACTORY_HPP
#define PLAYA_ILUKPRECONDITIONERFACTORY_HPP

#include "PlayaDefs.hpp"
#include "PlayaPreconditionerFactoryBase.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaILUFactorizableOp.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"

namespace Playa
{
  using namespace Teuchos;

  /**
   * 
   */
  template <class Scalar>
  class ILUKPreconditionerFactory
    : public PreconditionerFactoryBase<Scalar>
  {
  public:
    /** Construct with a parameter list */
    ILUKPreconditionerFactory(const ParameterList& params)
      : fillLevels_(1),
        overlapFill_(0),
        relaxationValue_(0.0),
        relativeThreshold_(1.0),
        absoluteThreshold_(0.0),
        leftOrRight_(Right)
    {
      LinearSolverBase<Scalar>::template setParameter<int>(params, &fillLevels_, 
                                                  "Graph Fill");

      LinearSolverBase<Scalar>::template setParameter<int>(params, &overlapFill_, 
                                                  "Overlap");

      LinearSolverBase<Scalar>::template setParameter<double>(params, &relaxationValue_, 
                                                     "Relaxation");

      LinearSolverBase<Scalar>::template setParameter<double>(params, &absoluteThreshold_, 
                                                     "Absolute Threshold");

      LinearSolverBase<Scalar>::template setParameter<double>(params, &relativeThreshold_, 
                                                     "Relative Threshold");

      bool isLeft = false;

      LinearSolverBase<Scalar>::template setParameter<bool>(params, &isLeft, "Left");

      if (isLeft) leftOrRight_ = Left;
      
    }


    /** virtual dtor */
    virtual ~ILUKPreconditionerFactory(){;}

    
    /** */
    virtual Preconditioner <Scalar>
    createPreconditioner(const LinearOperator<Scalar>& A) const 
    {
      /* In order for ILU factorization to work, the operator A must
       * implement the ILUFactorizableOp interface. We cast A's pointer
       * to a ILUFactorizableOp ptr. If the cast fails, throw a spoke. */
      
      const ILUFactorizableOp<Scalar>* fop 
        = dynamic_cast<const ILUFactorizableOp<Scalar>*>(A.ptr().get());

      TEST_FOR_EXCEPTION(fop==0, std::runtime_error,
                         "ILUKPreconditionerFactory attempted to "
                         "create an ILU preconditioner for an operator type "
                         "that does not implement the ILUFactorizableOp "
                         "interface. The op is " << A.description());

      
      /* Now we can delegate the construction of the ILU factors to 
      * the factorizable op. */
      Preconditioner<Scalar> P;
      fop->getILUKPreconditioner(fillLevels_,
                                 overlapFill_,
                                 relaxationValue_,
                                 relativeThreshold_,
                                 absoluteThreshold_,
                                 leftOrRight_,
                                 P);
      /* Return the preconditioner */
      return P;
    }

    /* Handleable boilerplate */
    GET_RCP(PreconditionerFactoryBase<Scalar>);
  private:

    int fillLevels_;
    int overlapFill_;
    Scalar relaxationValue_;
    Scalar relativeThreshold_;
    Scalar absoluteThreshold_;
    LeftOrRight leftOrRight_;
  };


}

#endif
