/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SUMOFBCS_H
#define SUNDANCE_SUMOFBCS_H

#include "SundanceDefs.hpp"
#include "SundanceSumOfIntegrals.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace FrameworkInterface;
  using std::string;

  namespace Internal
    {
      /** 
       * SumOfBCs represents a sum of essential
       * boundary conditions in integral form
       */
      class SumOfBCs : public SumOfIntegrals
        {
        public:
          /** Construct given an integral over a single domain */
          SumOfBCs(const RefCountPtr<CellFilterBase>& domain,
                         const Expr& expr,
                         const RefCountPtr<QuadratureFamilyBase>& quad);

          /** */
          virtual ~SumOfBCs(){;}

          /** Write a simple text description suitable 
           * for output to a terminal */
          virtual ostream& toText(ostream& os, bool paren) const ;

          /** Write in a form suitable for LaTeX formatting */
          virtual ostream& toLatex(ostream& os, bool paren) const ;

          /** Write in XML */
          virtual XMLObject toXML() const ;

          /** */
          virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

        private:
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
