/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ASSEMBLER_H
#define SUNDANCE_ASSEMBLER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceIntegrator.hpp"
#include "SundanceInserterBase.hpp"
#include "SundanceInserterFactoryBase.hpp"
#include "SundanceEvalManager.hpp"

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
    class Assembler : public TSFExtended::ObjectWithVerbosity<Assembler>
    {
    public:
      /** */
      Assembler(const Mesh& mesh, 
                const RefCountPtr<EquationSet>& eqn,
                const RefCountPtr<InserterFactoryBase>& inserterFactory,
                const VectorType<double>& vectorType,
                const VerbositySetting& verb = classVerbosity());
      
      /** */
      const RefCountPtr<DOFMapBase>& rowMap() const 
      {return inserter_->rowMap();}

      /** */
      const RefCountPtr<DOFMapBase>& colMap() const 
      {return inserter_->colMap();}

      /** */
      const RefCountPtr<Set<int> >& bcRows() const 
      {return inserter_->bcRows();}

      /** */
      void assemble(LinearOperator<double>& A,
                    Vector<double>& b) const ;

      /** */
      static int& workSetSize() 
      {static int rtn = defaultWorkSetSize(); return rtn;}
      
      /** */
      void getGraph(Array<Set<int> >& graph) const ;
      
    private:

      /** */
      bool isBCRow(int dof) const {return isBCRow_[dof-lowestRow_];}

      /** */
      static int defaultWorkSetSize() {return 100;}
      

      Mesh mesh_;

      RefCountPtr<EquationSet> eqn_;

      Array<RegionQuadCombo> rqc_;

      Array<int> isBCRqc_;

      RefCountPtr<InserterBase> inserter_;

      Array<RefCountPtr<Integrator> > integrator_;

      RefCountPtr<EvalManager> evalMgr_;

      Array<int> isBCRow_;

      int lowestRow_;

    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
