/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_BELOSSOLVER_HPP
#define PLAYA_BELOSSOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"
#include "PlayaPreconditionerFactory.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "PlayaBelosAdapter.hpp"
#include "Teuchos_Describable.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosTFQMRSolMgr.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Playa
{
  using namespace Teuchos;

  /**
   *
   */
  class BelosSolver : public LinearSolverBase<double>,
                      public Handleable<LinearSolverBase<double> >,
                      public Printable,
                      public Describable
  {
  public:
    /** */
    BelosSolver(const Teuchos::ParameterList& params);

    /** */
    virtual ~BelosSolver(){;}

    /** Set the preconditioning operator */
    void setUserPrec(const PreconditionerFactory<double>& pf) {pf_=pf;}

    /** \name Printable interface */
    //@{
    /** Write to a stream  */
    void print(std::ostream& os) const 
    {
      os << description() << std::endl;
    }
    //@}
    
    /** \name Describable interface */
    //@{
    /** Write a brief description */
    std::string description() const {return "BelosSolver";}
    //@}

    

    /** */
    virtual SolverState<double> solve(const LinearOperator<double>& op,
                                      const Vector<double>& rhs,
                                      Vector<double>& soln) const ;

    /** \name Handleable interface */
    //@{
    /** Return a ref count pointer to a newly created object */
    virtual RCP<LinearSolverBase<double> > getRcp() 
    {return rcp(this);}
    //@}


  protected:

  private:
    
    /** */
    PreconditionerFactory<double> pf_;
    /** */
    mutable RCP<Belos::SolverManager<double,Anasazi::SimpleMV, LinearOperator<double> > > solver_ ;
    /** */
    mutable bool hasSolver_;
  };
  
}

#endif
