/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_AMESOSSOLVER_HPP
#define PLAYA_AMESOSSOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "Teuchos_Describable.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
  using namespace Teuchos;

  /**
   *
   */
  class DenseLUSolver : public LinearSolverBase<double>,
                        public Playa::Handleable<LinearSolverBase<double> >,
                        public Printable,
                        public Describable
  {
  public:
    /** */
    DenseLUSolver();

    /** */
    virtual ~DenseLUSolver(){;}

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
    std::string description() const {return "DenseLUSolver";}
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
  };
  
}

#endif
