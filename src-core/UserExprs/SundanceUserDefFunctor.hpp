/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_USERDEFFUNCTOR_H
#define SUNDANCE_USERDEFFUNCTOR_H

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "Teuchos_Array.hpp"



namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

    /**
     * UserDefFunctor defines an interface for callbacks used to implement
     * user-defined nonlinear operators in the Sundance Expr system.
     */
    class UserDefFunctor
    {
    public:
      /** ctor */
      UserDefFunctor(const string& name) : name_(name) {;}

      /** */
      virtual ~UserDefFunctor(){;}

      /** */
      const string& name() const {return name_;}

      /** */
      virtual double eval0(const Array<double>& vars) const = 0 ;
                        

    private:
      string name_;
    };


}


#endif
