/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_UNARYFUNCTOR_H
#define SUNDANCE_UNARYFUNCTOR_H

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  namespace Internal
  {
    /**
     * 
     */
    class UnaryFunctor
    {
    public:
      /** ctor */
      UnaryFunctor(const string& name) : name_(name) {;}

      /** */
      const string& name() const {return name_;}

      /** */
      virtual void eval(const double* const x, int nx, double* f) const = 0 ;

      /** */
      virtual void eval(const double* const x, 
                        int nx, 
                        double* f, 
                        double* df) const = 0 ;

      /** Specify whether we should test for NAN or INFINITE results. */
      static bool& checkResults() {static bool rtn = false; return rtn;}
    private:
      string name_;
    };


  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
