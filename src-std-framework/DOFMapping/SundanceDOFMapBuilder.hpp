/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DOFMAPBUILDER_H
#define SUNDANCE_DOFMAPBUILDER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceBasisFamily.hpp"
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
    class DOFMapBuilder : public TSFExtended::ObjectWithVerbosity<DOFMapBase>
    {
    public:
      /** */
      DOFMapBuilder();
      /** */
      DOFMapBuilder(const Mesh& mesh, const RefCountPtr<EquationSet>& eqn);

      /** */
      const RefCountPtr<DOFMapBase>& rowMap() const {return rowMap_;}

      /** */
      const RefCountPtr<DOFMapBase>& colMap() const {return colMap_;}

      /** */
      const RefCountPtr<Array<int> >& isBCRow() const {return isBCRow_;}

      Array<BasisFamily> testBasisArray() const ;

      Array<BasisFamily> unkBasisArray() const ;

      const Mesh& mesh() const {return mesh_;}

    private:

      bool hasUnks() const ;

      bool unksAreHomogeneous() const ;

      bool testsAreHomogeneous() const ;

      bool unksAreOmnipresent() const ;

      bool testsAreOmnipresent() const ;

      bool regionIsMaximal(int r) const ;

      bool isSymmetric() const ;

      void markBCRows() ;

      const MPIComm& comm() const {return mesh().comm();}

      void init();

      Mesh mesh_;

      RefCountPtr<EquationSet> eqn_;

      RefCountPtr<DOFMapBase> rowMap_;

      RefCountPtr<DOFMapBase> colMap_;

      RefCountPtr<Array<int> > isBCRow_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
