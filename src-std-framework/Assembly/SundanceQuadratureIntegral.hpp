/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_QUADRATUREINTEGRAL_H
#define SUNDANCE_QUADRATUREINTEGRAL_H

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
                         int testDerivOrder,
                         const QuadratureFamily& quad);

      /** Construct a two-form to be computed by quadrature */
      QuadratureIntegral(int dim,
                         const CellType& cellType,
                         const BasisFamily& testBasis,
                         int testDerivOrder,
                         const BasisFamily& unkBasis,
                         int unkDerivOrder,
                         const QuadratureFamily& quad);
      
      /** */
      void transformTwoForm(const CellJacobianBatch& J, 
                            const Array<int>& alpha,
                            const Array<int>& beta,
                            const double* const coeff,
                            RefCountPtr<Array<double> >& A) const ;

      
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
