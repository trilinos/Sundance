/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_STDFWKEVALMEDIATOR_H
#define SUNDANCE_STDFWKEVALMEDIATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "SundanceDiscreteFunction.hpp"

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
     * StdFwkEvalMediator evaluates mesh-dependent functions in the
     * standard framework. A number of subtypes are supported: 
     * QuadratureEvalMediator, which does evaluation on quadrature points,
     * and NodalEvalMediator, which does evaluation at nodal points. 
     */
    class StdFwkEvalMediator : public AbstractEvalMediator,
                               public ObjectWithVerbosity<StdFwkEvalMediator>,
                               public TSFExtended::Printable
    {
    public:
      /** */
      StdFwkEvalMediator(const Mesh& mesh, int cellDim)
        : mesh_(mesh),
          cellDim_(cellDim),
          cellType_(NullCell),
          cellLID_(),
          cacheIsValid_(false)
      {;}

      /** */
      virtual ~StdFwkEvalMediator(){;}

      /** */
      void setCellBatch(const RefCountPtr<Array<int> >& cellLID) ;

      /** */
      virtual void setCellType(const CellType& cellType) 
      {cellType_=cellType; cacheIsValid() = false; jCacheIsValid_=false;}


    protected:
      const Mesh& mesh() const {return mesh_;}

      int cellDim() const {return cellDim_;}

      const CellType& cellType() const {return cellType_;}

      const RefCountPtr<Array<int> >& cellLID() const {return cellLID_;}

      bool& cacheIsValid() const {return cacheIsValid_;}


      /** */
      Map<const DiscreteFunction*, RefCountPtr<Array<double> > >& fCache() const {return fCache_;}
      /** */
      Map<const DiscreteFunction*, RefCountPtr<Array<double> > >& dfCache() const {return dfCache_;}
      /** */
      Map<const DiscreteFunction*, bool>& fCacheIsValid() const {return fCacheIsValid_;}
      /** */
      Map<const DiscreteFunction*, bool>& dfCacheIsValid() const {return dfCacheIsValid_;}
      
    private:
      Mesh mesh_;

      int cellDim_;

      CellType cellType_;

      RefCountPtr<Array<int> > cellLID_;

      mutable RefCountPtr<CellJacobianBatch> J_;

      mutable bool cacheIsValid_;

      mutable bool jCacheIsValid_;

      /** */
      mutable Map<const DiscreteFunction*, RefCountPtr<Array<double> > > fCache_;
      /** */
      mutable Map<const DiscreteFunction*, RefCountPtr<Array<double> > > dfCache_;

      /** */
      mutable Map<const DiscreteFunction*, bool> fCacheIsValid_;
      /** */
      mutable Map<const DiscreteFunction*, bool> dfCacheIsValid_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
