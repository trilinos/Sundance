/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_GROUPER_H
#define SUNDANCE_GROUPER_H

#include "SundanceDefs.hpp"
#include "SundanceSparsityPattern.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceEquationSet.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * Grouper
     */
    class GrouperBase
      : public TSFExtended::ObjectWithVerbosity<GrouperBase>
    {
    public:
      /** */
      GrouperBase() {;}

      /** */
      virtual ~GrouperBase(){;}

      /** */
      virtual void findGroups(const EquationSet& eqn,
                              const CellType& cellType,
                              int cellDim,
                              const QuadratureFamily& quad,
                              const RefCountPtr<SparsityPattern>& sparsity,
                              Array<IntegralGroup>& groups) const = 0 ;

    protected:
      void extractWeakForm(const EquationSet& eqn,
                           const MultipleDeriv& functionalDeriv,
                           BasisFamily& testBasis, 
                           BasisFamily& unkBasis,
                           MultiIndex& miTest, MultiIndex& miUnk,
                           int& testID, int& unkID, 
                           bool& isOneForm) const ;
                              
    };

  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
