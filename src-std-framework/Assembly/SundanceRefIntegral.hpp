/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_REFINTEGRAL_H
#define SUNDANCE_REFINTEGRAL_H

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
     * RefIntegral represents the integrals of a product 
     * of basis functions (or their derivatives) over a reference cell. 
     * This object can be created once and for all up front, and then 
     * reused to produce the integrals of constant-coefficient weak forms
     * over simplicial cells by means of a linear transformation.
     * 
     * An instance of this object can represent either a one form,
     * \f[
     * W_{i\gamma} = \int_{T_R} D_\gamma \psi_i 
     * \f]  
     * or a two form,
     * \f[
     * W_{(ij)(\gamma\delta)} = \int_{T_R} D_\gamma \psi_i D_\delta \phi_j.
     * \f]  
     * In the notation for the two form, we have grouped together the index
     * pairs \f$(ij)\f$ and \f$(\gamma\delta)\f$ to indicate that we 
     * will think of this 4-tensor as a 2D array, with \f$(ij)\f$ combinations
     * running down rows and \f$(\gamma\delta)\f$ running across columns. 
     * We will consider the node number combination 
     * \f$(ij)\f$ to be a single index \f$k\f$, and 
     * \f$(\gamma\delta)\f$ to be a single index \f$\epsilon\f$. Thus, 
     * the storage and use of one-forms and two-forms is essentially identical.
     * 
     * This object's job in life is to be multiplied with a transformation
     * matrix \f$T_{\epsilon c}\f$to produce an array of local element matrices
     * \f$A_{kc}\f$,
     * \f[
     * A_{kc} = W_{k\epsilon} T_{\epsilon c}
     * \f]
     * The index \f$c\f$ is over cells (cells are processed in batches).
     * 
     * Physical storage is as a 1D vector stored in colum-major order. 
     */
    class RefIntegral : public ElementIntegral
    {
    public:
      /** Construct a reference zero-form */
      RefIntegral(int dim, 
                  const CellType& cellType);

      /** Construct a reference one-form */
      RefIntegral(int dim, 
                  const CellType& cellType,
                  const BasisFamily& testBasis,
                  const Array<int>& alpha,
                  int testDerivOrder);

      /** Construct a reference two-form */
      RefIntegral(int dim,
                  const CellType& cellType,
                  const BasisFamily& testBasis,
                  const Array<int>& alpha,
                  int testDerivOrder,
                  const BasisFamily& unkBasis,
                  const Array<int>& beta,
                  int unkDerivOrder);
      
      /** */
      void print(ostream& os) const ;

      /** */
      void transform(const CellJacobianBatch& J, 
                     const Array<double>& coeff,
                     RefCountPtr<Array<double> >& A) const 
      {
        if (order()==2) transformTwoForm(J, coeff, A);
        else if (order()==1) transformOneForm(J, coeff, A);
        else transformZeroForm(J, coeff, A);
      }

      /** */
      void transformTwoForm(const CellJacobianBatch& J, 
                            const Array<double>& coeff,
                            RefCountPtr<Array<double> >& A) const ;

      /** */
      void transformOneForm(const CellJacobianBatch& J, 
                            const Array<double>& coeff,
                            RefCountPtr<Array<double> >& A) const ;

      /** */
      void transformZeroForm(const CellJacobianBatch& J, 
                            const Array<double>& coeff,
                            RefCountPtr<Array<double> >& A) const ;
      /** */
      inline double& value(int testDerivDir, int testNode,
                           int unkDerivDir, int unkNode)
      {return W_[unkNode + nNodesUnk()*testNode 
                 + nNodes()*(unkDerivDir + nRefDerivUnk()*testDerivDir)];}

      

      /** */
      inline const double& value(int testDerivDir, int testNode,
                                 int unkDerivDir, int unkNode) const 
      {
        return W_[unkNode + nNodesUnk()*testNode 
                  + nNodes()*(unkDerivDir + nRefDerivUnk()*testDerivDir)];
      }
      
      /** */
      inline double& value(int testDerivDir, int testNode)
      {return W_[nNodesTest()*testDerivDir + testNode];}


      /** */
      inline const double& value(int testDerivDir, int testNode) const 
      {return W_[nNodesTest()*testDerivDir + testNode];}

      static double& totalFlops() {static double rtn = 0; return rtn;}

    protected:

      static void addFlops(const double& flops) {totalFlops() += flops;}

      /** */
      void createTwoFormTransformationMatrix(const CellJacobianBatch& J,  
                                             const Array<int>& alpha,
                                             const Array<int>& beta,
                                             const Array<double>& coeff) const;
      /** */
      void createOneFormTransformationMatrix(const CellJacobianBatch& J,  
                                             const Array<int>& alpha,
                                             const Array<double>& coeff) const;
      
    private:

      Array<double> W_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
