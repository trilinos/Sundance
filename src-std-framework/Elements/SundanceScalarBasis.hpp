/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SCALARBASIS_H
#define SUNDANCE_SCALARBASIS_H

#include "SundanceDefs.hpp"
#include "SundanceBasisFamilyBase.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace Teuchos;
  using namespace TSFExtended;
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  
  namespace Internal
  {
    /** 
     * Base class for scalar-valued basis families
     */
    class ScalarBasis : public BasisFamilyBase
    {
    public:
      /** */
      ScalarBasis(){;}

      /** */
      virtual ~ScalarBasis(){;}

      /** Return the spatial dimension of the basis */
      virtual int dim() const {return 1;}
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
