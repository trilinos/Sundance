/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_NONLINEAREXPR_H
#define SUNDANCE_NONLINEAREXPR_H

#include "SundanceDefs.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceSet.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using std::string;

  namespace Internal
    {
      /** 
       *
       */
      class NonlinearExpr
      {
      public:
        /** */
        NonlinearExpr(){;}
        
      protected:
        
        Set<MultiSet<int> > 
        findChildFuncIDSet(const Set<MultiSet<int> >& activeFuncIDs,
                           const Set<int>& allFuncIDs) const ;
      };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
