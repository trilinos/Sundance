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
      RefBasisProductIntegral(int testDerivOrder,
                              int nRefDerivTest, 
                              int nNodesTest,
                              const Array<int>& alpha,
                              const Array<double>& coeff)
        : W_(nRefDerivTest*nNodesTest),
          testDerivOrder_(testDerivOrder), 
          nRefDerivTest_(nRefDerivTest),
          nNodesTest_(nNodesTest),
          unkDerivOrder_(-1), 
          nRefDerivUnk_(-1),
          nNodesUnk_(-1),
          alpha_(alpha),
          beta_(),
          coeff_(coeff),
          isTwoForm_(false)
      {;}

      /** Construct a reference two-form */
      RefBasisProductIntegral(int testDerivOrder,
                              int nRefDerivTest, 
                              int nNodesTest,
                              int unkDerivOrder,
                              int nRefDerivUnk, 
                              int nNodesUnk,
                              const Array<int>& alpha,
                              const Array<int>& beta,
                              const Array<double>& coeff)
        : W_(nRefDerivTest*nNodesTest*nRefDerivUnk*nNodesUnk),
          testDerivOrder_(testDerivOrder), 
          nRefDerivTest_(nRefDerivTest),
          nNodesTest_(nNodesTest), 
          unkDerivOrder_(unkDerivOrder), 
          nRefDerivUnk_(nRefDerivUnk),
          nNodesUnk_(nNodesUnk), 
          alpha_(alpha),
          beta_(beta),
          coeff_(coeff),
          isTwoForm_(true)
      {;}
      
      /** */
      void transformToPhysicalCoords(const CellJacobianBatch& J, 
                                     RefCountPtr<Array<double> >& A) const ;

          
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

      int testDerivOrder_;

      int nRefDerivTest_;

      int nNodesTest_;

      int unkDerivOrder_;

      int nRefDerivUnk_;

      int nNodesUnk_;

      Array<int> alpha_;

      Array<int> beta_;

      Array<double> coeff_;

      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
