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
     * 
     */
    class UserDefFunctor
    {
    public:
      /** ctor */
      UserDefFunctor(const string& name) : name_(name) {;}

      /** */
      const string& name() const {return name_;}

      /** */
      virtual double eval(const Array<double>& vars) const = 0 ;
                        

    private:
      string name_;
    };


}


#endif
