/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_INTEGRATOR_H
#define SUNDANCE_INTEGRATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceLocalMatrixContainer.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "TSFObjectWithVerbosity.hpp"

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
     * Integrator
     */
    class Integrator
      : public TSFExtended::ObjectWithVerbosity<Integrator>
    {
    public:
      /** */
      Integrator(const Mesh& mesh, 
                 const Expr& expr,
                 const DerivSet& nonzeroDerivs,
                 const RegionQuadCombo& rqc,
                 const RefCountPtr<EvalManager>& evalMgr);
      
      /** */
      void integrate(const RefCountPtr<Array<int> >& workSet,
                     RefCountPtr<LocalMatrixContainer>& localMat) const ;

    protected:
      /** */
      int cellDim() const {return cellDim_;}

      /** */
      const CellType& cellType() const {return cellType_;}

      /** */
      const Mesh& mesh() const {return mesh_;}

      /** */
      const EvaluatableExpr *const expr() const {return expr_;}

      /** */
      const RefCountPtr<SparsitySuperset>& sparsity() const {return sparsity_;}

    private:
      /** */
      int cellDim_;

      /** */
      CellType cellType_;

      /** */
      Mesh mesh_;

      /** */
      const EvaluatableExpr *const expr_;

      /** */
      mutable RefCountPtr<EvalManager> evalMgr_;

      /** */
      DerivSet nonzeroDerivs_;

      /** */
      RegionQuadCombo rqc_;

      /** */
      mutable RefCountPtr<StdFwkEvalMediator> mediator_;

      /** */
      RefCountPtr<SparsitySuperset> sparsity_;

    };

  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
