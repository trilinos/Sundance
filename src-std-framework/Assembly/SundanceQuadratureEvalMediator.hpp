/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_QUADRATUREEVALMEDIATOR_H
#define SUNDANCE_QUADRATUREEVALMEDIATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceOrderedTuple.hpp"

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
    class QuadratureEvalMediator : public StdFwkEvalMediator
    {
    public:
      /** 
       * 
       */
      QuadratureEvalMediator(const Mesh& mesh, 
                             int cellDim,
                             const QuadratureFamily& quad);

      /** */
      virtual ~QuadratureEvalMediator(){;}

      
      /** Evaluate the given coordinate expression, putting
       * its numerical values in the given EvalVector. */
      virtual void evalCoordExpr(const CoordExpr* expr,
                                 RefCountPtr<EvalVector>& vec) const ;
      
      /** Evaluate the given discrete function, putting
       * its numerical values in the given EvalVector. */
      virtual void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                           const Array<MultiIndex>& mi,
                                           Array<RefCountPtr<EvalVector> >& vec) const ;

      /** Evaluate the given cell diameter expression, putting
       * its numerical values in the given EvalVector. */
      virtual void evalCellDiameterExpr(const CellDiameterExpr* expr,
                                        RefCountPtr<EvalVector>& vec) const ;

      /** */
      virtual void setCellType(const CellType& cellType) ;

      /** */
      virtual void print(ostream& os) const ;

      /** */
      RefCountPtr<Array<Array<Array<double> > > > 
      getRefBasisVals(const BasisFamily& basis, 
                      int diffOrder) const ;

      /** */
      const Array<double>& quadWgts() const 
      {return *(refQuadWeights_.get(cellType()));}

      static double& totalFlops() {static double rtn = 0; return rtn;}

      

      static void addFlops(const double& flops) {totalFlops() += flops;}

    private:


      /** */
      void fillFunctionCache(const DiscreteFunction* f,
                             const MultiIndex& mi) const ;

     
      /** */
      void computePhysQuadPts() const ;

      /** */
      QuadratureFamily quad_;

      /** */
      Map<CellType, RefCountPtr<Array<Point> > > refQuadPts_;

      /** */
      Map<CellType, RefCountPtr<Array<double> > > refQuadWeights_;

      /** */
      mutable Array<Point> physQuadPts_;

      /** */
      mutable Array<Map<OrderedPair<BasisFamily, CellType>, RefCountPtr<Array<Array<Array<double> > > > > > refBasisVals_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
