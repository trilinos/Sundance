/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_REFBASISPRODUCTINTEGRAL_H
#define SUNDANCE_REFBASISPRODUCTINTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceJacobianBatch.hpp"
#include "Teuchos_Array.hpp"

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
     * RefBasisProductIntegral represents the integrals of a product 
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
    class RefBasisProductIntegral
    {
    public:
      /** Construct a reference one-form */
      RefBasisProductIntegral(int dim, 
                              const CellType& cellType,
                              const BasisFamily& testBasis,
                              int testOrder,
                              const QuadratureFamily& quad);

      /** Construct a reference two-form */
      RefBasisProductIntegral(int dim,
                              const CellType& cellType,
                              const BasisFamily& testBasis,
                              int testOrder,
                              const BasisFamily& unkBasis,
                              int unkOrder,
                              const QuadratureFamily& quad);
        

      /** */
      bool isTwoForm() const {return isTwoForm_;}
      
      /** */
      void transformToPhysicalCoords(const CellJacobianBatch& J, 
                                     const Array<int>& alpha,
                                     const Array<int>& beta,
                                     const Array<double>& coeff,
                                     RefCountPtr<Array<double> >& A) const ;

      
      /** */
      int nRefDerivTest() const {return nRefDerivTest_;}
      
      /** */
      int nRefDerivUnk() const {return nRefDerivUnk_;}


      /** */
      inline double& value(int testDerivDir, int testNode,
                           int unkDerivDir, int unkNode)
      {return W_[testNode + nNodesTest_*unkNode 
                    + nNodes_*(unkDerivDir + nRefDerivsUnk_*testDerivDir)];}

      

      /** */
      inline const double& value(int testDerivDir, int testNode,
                                 int unkDerivDir, int unkNode) const 
      {
        return W_[testNode + nNodesTest_*unkNode 
                     + nNodes_*(unkDerivDir + nRefDerivsUnk_*testDerivDir)];
      }
      
      /** */
      inline double& value(int testDerivDir, int testNode)
      {return W_[nNodesTest_*testDerivDir + testNode];}


      /** */
      inline const double& value(int testDerivDir, int testNode) const 
      {return W_[nNodesTest_*testDerivDir + testNode];}

    protected:

      /** */
      static Array<int>& G() {static Array<int> rtn; return rtn;}

      /** */
      static Array<int>& T() {static Array<int> rtn; return rtn;}

    private:

      Array<double> W_;

      int dim_;

      int testOrder_;

      int nRefDerivTest_;

      int nNodesTest_;

      int unkOrder_;

      int nRefDerivUnk_;

      int nNodesUnk_;

      bool isTwoForm_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
