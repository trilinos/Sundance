/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EXPRBASE_H
#define SUNDANCE_EXPRBASE_H


#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  namespace Internal
    {
      /** */
      class ExprBase
        {
        public:
          /** empty ctor */
          ExprBase();

          /** virtual destructor */
          virtual ~ExprBase() {;}

          /** Write a simple text description suitable 
           * for output to a terminal */
          virtual ostream& toText(ostream& os, bool paren) const = 0 ;

          /** Write in a form suitable for LaTeX formatting */
          virtual ostream& toLatex(ostream& os, bool paren) const = 0 ;

          /** Append to the set of unknown func IDs present in this expression.
           * Base class does nothing */
          virtual void accumulateUnkSet(Set<int>& unkIDs) const {;}

          /** Append to the set of test func IDs present in this expression.
           * Base class does nothing */
          virtual void accumulateTestSet(Set<int>& testIDs) const {;}

          /** */
          string toString() const ;

          /** Write in XML */
          virtual XMLObject toXML() const = 0 ;

          /** Return a descriptive name for the expression subtype */
          virtual string typeName() const ;

          /** */
          virtual RefCountPtr<ExprBase> getRcp() = 0 ;

        protected:
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
