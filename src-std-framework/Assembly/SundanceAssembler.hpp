/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ASSEMBLER_H
#define SUNDANCE_ASSEMBLER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceGrouperBase.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "TSFLoadableVector.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFVector.hpp"
#include "TSFVectorType.hpp"

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
     * 
     */
    class Assembler : public TSFExtended::ObjectWithVerbosity<Assembler>
    {
    public:
      /** */
      Assembler(const Mesh& mesh, 
                const RefCountPtr<EquationSet>& eqn,
                const VectorType<double>& vectorType,
                const VerbositySetting& verb = classVerbosity());
      
      /** */
      const RefCountPtr<DOFMapBase>& rowMap() const 
      {return rowMap_;}

      /** */
      const RefCountPtr<DOFMapBase>& colMap() const 
      {return colMap_;}

      /** */
      const RefCountPtr<Set<int> >& bcRows() {return bcRows_;}

      /** */
      void assemble(TSFExtended::LinearOperator<double>& A,
                    TSFExtended::Vector<double>& b) const ;

      /** */
      static int& workSetSize() 
      {static int rtn = defaultWorkSetSize(); return rtn;}
      
      /** */
      void getGraph(Array<Set<int> >& graph) const ;
      
    private:

      /** */
      void insertLocalMatrixValues(int cellDim, const Array<int>& workSet,
                                   bool isBCRqc, 
                                   int nTestNodes, int nUnkNodes,
                                   const Array<int>& testID,
                                   const Array<int>& unkID, 
                                   const Array<double>& localValues,
                                   TSFExtended::LoadableMatrix<double>* mat) const ;

      /** */
      void insertLocalVectorValues(int cellDim, const Array<int>& workSet,
                                   bool isBCRqc, 
                                   int nTestNodes, 
                                   const Array<int>& testID,
                                   const Array<double>& localValues,
                                   TSFExtended::LoadableVector<double>* vec) const ;

      /** */
      bool isBCRow(int dof) const {return isBCRow_[dof-lowestRow_];}

      /** */
      static int defaultWorkSetSize() {return 100;}
      

      Mesh mesh_;

      RefCountPtr<EquationSet> eqn_;

      RefCountPtr<DOFMapBase> rowMap_;

      RefCountPtr<DOFMapBase> colMap_;

      RefCountPtr<Set<int> > bcRows_;

      Array<RegionQuadCombo> rqc_;

      Array<int> isBCRqc_;

      Array<Array<IntegralGroup> > groups_;

      Array<RefCountPtr<StdFwkEvalMediator> > mediators_;

      Array<const EvaluatableExpr*> evalExprs_;

      RefCountPtr<EvalManager> evalMgr_;

      Array<int> isBCRow_;

      int lowestRow_;

      VectorType<double> vecType_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
