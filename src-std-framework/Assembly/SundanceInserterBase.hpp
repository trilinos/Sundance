/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_INSERTERBASE_H
#define SUNDANCE_INSERTERBASE_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceLocalMatrixContainer.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "TSFVectorType.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFVector.hpp"

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
     * InserterBase provides an interface for insertion of local matrix and
     * vector values in a global matrix and vector
     */
    class InserterBase
      : public TSFExtended::ObjectWithVerbosity<InserterBase>
    {
    public:
      /** */
      InserterBase(const RefCountPtr<DOFMapBase>& rowMap,
                   const RefCountPtr<DOFMapBase>& colMap,
                   const RefCountPtr<Set<int> >& bcRows) 
        : rowMap_(rowMap),
          colMap_(colMap),
          bcRows_(bcRows)
      {;}

      /** virtual dtor */
      virtual ~InserterBase() {;}

      /** */
      virtual void insert(int cellDim, const RefCountPtr<Array<int> >& cells,
                          bool isBC,
                          const LocalMatrixContainer* localMat,
                          LinearOperator<double>& A,
                          Vector<double>& b) const = 0 ;

      /** */
      virtual void configureMat(LinearOperator<double>& A) const {;}

      

      /** */
      const RefCountPtr<DOFMapBase>& rowMap() const {return rowMap_;}
      /** */
      const RefCountPtr<DOFMapBase>& colMap() const {return colMap_;}
      /** */
      const RefCountPtr<Set<int> >& bcRows() const {return bcRows_;}
    private:
      RefCountPtr<DOFMapBase> rowMap_;
      RefCountPtr<DOFMapBase> colMap_;
      RefCountPtr<Set<int> > bcRows_;
    };

  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
