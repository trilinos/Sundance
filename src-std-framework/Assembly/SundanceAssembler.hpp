/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ASSEMBLER_H
#define SUNDANCE_ASSEMBLER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceWeakFormBatch.hpp"
#include "TSFObjectWithVerbosity.hpp"

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
    class Assembler : public TSFExtended::ObjectWithVerbosity<Assembler>,
                      public TSFExtended::Printable
    {
    public:
      /** */
      Assembler(const Mesh& mesh, const RefCountPtr<EquationSet>& eqn);

      /** */
      const RefCountPtr<DOFMapBase>& rowMap() const {return rowMap_;}

      /** */
      const RefCountPtr<DOFMapBase>& colMap() const {return colMap_;}

      /** */
      const RefCountPtr<Set<int> >& bcRows() const {return bcRows_;}

      /** */
      void print(ostream& os) const ;

    private:
      
      void addToWeakFormBatch(const DerivSet& derivs);

      void getGraph(Array<Set<int> >& graph);

      Mesh mesh_;

      RefCountPtr<EquationSet> eqn_;

      RefCountPtr<DOFMapBase> rowMap_;

      RefCountPtr<DOFMapBase> colMap_;

      RefCountPtr<Set<int> > bcRows_;

      Array<RegionQuadCombo> rqc_;

      Array<int> isBCRqc_;

      Array<Expr> rqcExprs_;

      Array<DerivSet> rqcDerivSet_;

      Array<Array<RefCountPtr<WeakFormBatch> > > weakForms_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
