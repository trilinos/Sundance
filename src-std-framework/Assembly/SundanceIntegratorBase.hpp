/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_INTEGRATORBASE_H
#define SUNDANCE_INTEGRATORBASE_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceLocalMatrixBatch.hpp"
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
     * IntegratorBase provides an interface and common implementation
     * for objects that carry out integrations on elements. 
     */
    class IntegratorBase 
      : public TSFExtended::ObjectWithVerbosity<IntegratorBase>
    {
    public:
      /** */
      IntegratorBase(const Mesh& mesh, 
                     const Expr& expr,
                     const DerivSet& nonzeroDerivs,
                     const RegionQuadCombo& rqc,
                     const RefCountPtr<EvalManager>& evalMgr);
      
      /** */
      void integrate(const RefCountPtr<Array<int> >& workSet,
                     RefCountPtr<LocalMatrixBatch>& localMat) const ;

    protected:
      /** */
      virtual RefCountPtr<StdFwkEvalMediator> 
      createEvalMediator(const Mesh& mesh, const RegionQuadCombo& rqc) const = 0 ;

      /** */
      virtual void innerIntegrate(const RefCountPtr<Array<int> >& workSet,
                                  RefCountPtr<LocalMatrixBatch>& localMat) const = 0 ;

      /** */
      virtual void init() = 0 ;

      /** */
      int cellDim() const {return cellDim_;}

      /** */
      const CellType& cellType() const {return cellType_;}

      /** */
      const Mesh& mesh() const {return mesh_;}

      /** */
      const EvaluatableExpr *const expr() const {return expr_;}

      /** */
      void evaluate(RefCountPtr<EvalVectorArray>& results) const ;

      /** */
      void getJacobians(const Array<int>& workSet,
                        CellJacobianBatch& J) const ;

      /** */
      const RefCountPtr<SparsityPattern>& sparsity() const {return sparsity_;}

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
      RefCountPtr<SparsityPattern> sparsity_;

      /** */
      mutable bool needsInit_;
    };

  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
