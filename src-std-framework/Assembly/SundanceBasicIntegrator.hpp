/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BASICINTEGRATOR_H
#define SUNDANCE_BASICINTEGRATOR_H

#include "SundanceDefs.hpp"
#include "SundanceIntegratorBase.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace TSFExtended;

  namespace Internal
  {
    using namespace Teuchos;

    /** 
     *
     */
    class BasicIntegrator : public IntegratorBase
    {
    public:
      /** */
      BasicIntegrator(const Mesh& mesh, 
                      const Expr& expr,
                      const DerivSet& nonzerDerivs,
                      const RegionQuadCombo& rqc,
                      const RefCountPtr<EvalManager>& evalMgr);
      
    protected:

      /** */
      virtual RefCountPtr<StdFwkEvalMediator> 
      createEvalMediator(const Mesh& mesh, const RegionQuadCombo& rqc) const ;

      /** */
      virtual void innerIntegrate(const RefCountPtr<Array<int> >& workSet,
                                  RefCountPtr<LocalMatrixContainer>& localMat) const ;

      /** */
      virtual void init() ;

      
    private:
    };

  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
