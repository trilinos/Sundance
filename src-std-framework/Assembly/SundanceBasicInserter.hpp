/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BASICINSERTER_H
#define SUNDANCE_BASICINSERTER_H

#include "SundanceDefs.hpp"
#include "SundanceInserterBase.hpp"

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
    class BasicInserter : public InserterBase
    {
    public:
      /** */
      BasicInserter(const RefCountPtr<DOFMapBase>& rowMap,
                   const RefCountPtr<DOFMapBase>& colMap,
                   const RefCountPtr<Set<int> >& bcRows) 
        : InserterBase(rowMap, colMap, bcRows)
      {;}

      /** virtual dtor */
      virtual ~BasicInserter() {;}

      /** */
      virtual void insert(int cellDim, const RefCountPtr<Array<int> >& cells,
                          bool isBC,
                          const LocalMatrixContainer* localMat,
                          LinearOperator<double>& A,
                          Vector<double>& b) const ;

      /** */
      virtual void configureMat(LinearOperator<double>& A) const ;
     
    };

  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
