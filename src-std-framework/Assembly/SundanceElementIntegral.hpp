/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ELEMENTINTEGRAL_H
#define SUNDANCE_ELEMENTINTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
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
     * ElementIntegral encapsulates the common data needed for the
     * integration of groups of related one-forms and two-forms.
     */
    class ElementIntegral 
      : public TSFExtended::ObjectWithVerbosity<ElementIntegral>,
        public TSFExtended::Printable
    {
    public:
      /** Construct a one-form */
      ElementIntegral(int dim, 
                      const CellType& cellType,
                      const BasisFamily& testBasis,
                      const Array<int>& alpha,
                      int testDerivOrder);

      /** Construct a two-form */
      ElementIntegral(int dim,
                      const CellType& cellType,
                      const BasisFamily& testBasis,
                      const Array<int>& alpha,
                      int testDerivOrder,
                      const BasisFamily& unkBasis,
                      const Array<int>& beta,
                      int unkDerivOrder);

      /** Indicate whether this element integral is a two form */
      bool isTwoForm() const {return isTwoForm_;}
      
      /** Return the number of nodes associated with the test function */
      int nNodesTest() const {return nNodesTest_;}
      
      /** Return the number of nodes associated with the test function */
      int nNodesUnk() const {return nNodesUnk_;}

      /** Return the total number of elements in this local stiffness
       * matrix */
      int nNodes() const {return nNodes_;}


    protected:

      /** The dimension of the cell being integrated */
      int dim() const {return dim_;}
      
      /** Number of test function derivatives wrt reference coordinates that
       * are needed to evaluate this integral. Will always be equal to 
       * ipow(element dimension, differentiation order). */
      int nRefDerivTest() const {return nRefDerivTest_;}
      
      /** Number of unknown function derivatives wrt reference coordinates that
       * are needed to evaluate this integral. Will always be equal to 
       * ipow(element dimension, differentiation order). */
      int nRefDerivUnk() const {return nRefDerivUnk_;}

      /** The order to which the test function is differentiated in this
       * integral. */
      int testDerivOrder() const {return testDerivOrder_;}

      /** The order to which the unknown function is differentiated in this
       * integral. */
      int unkDerivOrder() const {return unkDerivOrder_;}

      /** */
      const Array<int>& alpha() const {return alpha_;}

      /** */
      const Array<int>& beta() const {return beta_;}

      /** Workspace for element transformations */
      static Array<double>& G() {static Array<double> rtn; return rtn;}

      /** return base to the given power */
      static int ipow(int base, int power);

      /** The value below which chop() sets numbers to zero */
      static double chopVal() {static double rtn=1.0e-14; return rtn;}

      /** Chop a number to zero if it is smaller in magnitude than
       * the value chopVal() */
      static double chop(const double& x) 
      {
        if (::fabs(x) > chopVal()) return x;
        else return 0.0;
      }

    private:

      int dim_;

      int testDerivOrder_;

      int nRefDerivTest_;

      int nNodesTest_;

      int unkDerivOrder_;

      int nRefDerivUnk_;

      int nNodesUnk_;

      int nNodes_;

      bool isTwoForm_;

      Array<int> alpha_;

      Array<int> beta_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
