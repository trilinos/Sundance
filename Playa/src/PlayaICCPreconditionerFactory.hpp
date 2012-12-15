/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_ICCPRECONDITIONERFACTORY_HPP
#define PLAYA_ICCPRECONDITIONERFACTORY_HPP

#include "PlayaDefs.hpp"
#include "PlayaPreconditionerFactoryBase.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaICCFactorizableOp.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"

namespace Playa
{
  using namespace Teuchos;

  /**
   * 
   */
  template <class Scalar>
  class ICCPreconditionerFactory
    : public PreconditionerFactoryBase<Scalar>
  {
  public:
    /** Magnitude type */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    /** Construct with a parameter list */
    //const ParameterList& params inside
    ICCPreconditionerFactory():
    	fillLevels_(1),
        overlapFill_(0),
	dropTolerance_(1.0e-4),
        relaxationValue_(0.0),
        relativeThreshold_(1.0),
        absoluteThreshold_(0.0){;}
        
    ICCPreconditionerFactory(int fillLevels, 
			     int overlapFill, 
			     ScalarMag dropTolerance,
			     ScalarMag relaxationValue, 
			     ScalarMag relativeThreshold, 
			     ScalarMag absoluteThreshold)
    {
    	fillLevels_=fillLevels;
    	overlapFill_=overlapFill;
	dropTolerance_ = dropTolerance;
    	relaxationValue_=relaxationValue;
    	relativeThreshold_=relativeThreshold;
    	absoluteThreshold_=absoluteThreshold;
    }

    /** virtual dtor */
    virtual ~ICCPreconditionerFactory(){;}

    
    /** */
    virtual Preconditioner <Scalar>
    createPreconditioner(const LinearOperator<Scalar>& A) const 
    {
      /* In order for ICC factorization to work, the operator A must
       * implement the ICCFactorizableOp interface. We cast A's pointer
       * to a ICCFactorizableOp ptr. If the cast fails, throw a spoke. */
      
      const ICCFactorizableOp<Scalar>* fop 
        = dynamic_cast<const ICCFactorizableOp<Scalar>*>(A.ptr().get());

      TEUCHOS_TEST_FOR_EXCEPTION(fop==0, std::runtime_error,
                         "ICCPreconditionerFactory attempted to "
                         "create an ICC preconditioner for an operator type "
                         "that does not implement the ICCFactorizableOp "
                         "interface. The op is " << A.description());

      
      /* Now we can delegate the construction of the ICC factors to 
      * the factorizable op. */
      Preconditioner<Scalar> P;
      fop->getICCPreconditioner(fillLevels_,
				overlapFill_,
				dropTolerance_,
				relaxationValue_,
				relativeThreshold_,
				absoluteThreshold_,
				P);
      /* Return the preconditioner */
      return P;
    }

    /* Handleable boilerplate */
    GET_RCP(PreconditionerFactoryBase<Scalar>);
  private:

    int fillLevels_;
    int overlapFill_;
    ScalarMag dropTolerance_;
    ScalarMag relaxationValue_;
    ScalarMag relativeThreshold_;
    ScalarMag absoluteThreshold_;
  };


}

#endif
