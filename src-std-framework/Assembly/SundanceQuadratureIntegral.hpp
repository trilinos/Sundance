/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_QUADRATUREINTEGRAL_H
#define SUNDANCE_QUADRATUREINTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceElementIntegral.hpp"

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
     *
     */
    class QuadratureIntegral 
      : public ElementIntegral
    {
    public:
      /** Construct a one form to be computed by quadrature */
      QuadratureIntegral(int dim, 
                         const CellType& cellType,
                         const BasisFamily& testBasis,
                         const Array<int>& alpha,
                         int testDerivOrder,
                         const QuadratureFamily& quad);

      /** Construct a two-form to be computed by quadrature */
      QuadratureIntegral(int dim,
                         const CellType& cellType,
                         const BasisFamily& testBasis,
                         const Array<int>& alpha,
                         int testDerivOrder,
                         const BasisFamily& unkBasis,
                         const Array<int>& beta,
                         int unkDerivOrder,
                         const QuadratureFamily& quad);

      /** */
      void transform(const CellJacobianBatch& J, 
                     const double* const coeff,
                     RefCountPtr<Array<double> >& A) const 
      {
        if (isTwoForm()) transformTwoForm(J, coeff, A);
        else transformOneForm(J, coeff, A);
      }
      
      /** */
      void transformTwoForm(const CellJacobianBatch& J, 
                            const double* const coeff,
                            RefCountPtr<Array<double> >& A) const ;
      
      /** */
      void transformOneForm(const CellJacobianBatch& J, 
                            const double* const coeff,
                            RefCountPtr<Array<double> >& A) const ;

      /** */
      void print(ostream& os) const ;
      


      /** */
      int nQuad() const {return nQuad_;}

      static double& totalFlops() {static double rtn = 0; return rtn;}

    private:

      static void addFlops(const double& flops) {totalFlops() += flops;}

      /** Do the integration by summing reference quantities over quadrature
       * points and then transforming the sum to physical quantities.  */
      void transformSummingFirst(int nCells,
                                 const double* const coeff,
                                 RefCountPtr<Array<double> >& A) const ;

      /** Do the integration by transforming to physical coordinates 
       * at each quadrature point, and then summing */
      void transformSummingLast(int nCells,
                                const double* const coeff,
                                RefCountPtr<Array<double> >& A) const ;

      /** Determine whether to do this batch of integrals using the
       * sum-first method or the sum-last method */
      bool useSumFirstMethod() const {return useSumFirstMethod_;}
      
      /** */
      inline double& wValue(int q, int testDerivDir, int testNode,
                           int unkDerivDir, int unkNode)
      {return W_[unkNode
                  + nNodesUnk()
                  *(testNode + nNodesTest()
                    *(unkDerivDir + nRefDerivUnk()
                      *(testDerivDir + nRefDerivTest()*q)))];}

      

      /** */
      inline const double& wValue(int q, 
                                 int testDerivDir, int testNode,
                                 int unkDerivDir, int unkNode) const 
      {
        return W_[unkNode
                  + nNodesUnk()
                  *(testNode + nNodesTest()
                    *(unkDerivDir + nRefDerivUnk()
                      *(testDerivDir + nRefDerivTest()*q)))];
      }
      
      /** */
      inline double& wValue(int q, int testDerivDir, int testNode)
      {return W_[testNode + nNodesTest()*(testDerivDir + nRefDerivTest()*q)];}


      /** */
      inline const double& wValue(int q, int testDerivDir, int testNode) const 
      {return W_[testNode + nNodesTest()*(testDerivDir + nRefDerivTest()*q)];}

      /* */
      void createTwoFormTransformationMatrix(const CellJacobianBatch& J,  
                                             const Array<int>& alpha,
                                             const Array<int>& beta) const ;

      /* */
      void createOneFormTransformationMatrix(const CellJacobianBatch& J,  
                                             const Array<int>& alpha) const ;
      /* */
      Array<double> W_;

      /* */
      int nQuad_;

      /* */
      bool useSumFirstMethod_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
