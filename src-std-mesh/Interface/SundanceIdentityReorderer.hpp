/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_IDENTITYREORDERER_H
#define SUNDANCE_IDENTITYREORDERER_H


#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellReordererImplemBase.hpp"
#include "SundanceCellReordererBase.hpp"
#include "TSFHandleable.hpp"

namespace Sundance
{
  using namespace TSFExtended;
  namespace Internal
  {
    using namespace Teuchos;


    /**
     * The identity reorderer walks through cells in whatever
     * order they are numbered the the mesh. 
     */
    class IdentityReordererImplem : public CellReordererImplemBase
    {
    public:
      /** */
      IdentityReordererImplem(const MeshBase* mesh); 
      
      /** */
      virtual ~IdentityReordererImplem(){;}
    
      /** */
      virtual int advance(int currentLID) const {return currentLID+1;}
    };
  }

  /** */
  class IdentityReorderer 
    : public Internal::GenericCellReordererFactory<Internal::IdentityReordererImplem>
  {
  public:
    IdentityReorderer(){;}

    virtual ~IdentityReorderer(){;}

    GET_RCP(Internal::CellReordererFactoryBase);
  };
}

#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
