/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_INSERTERFACTORYBASE_H
#define SUNDANCE_INSERTERFACTORYBASE_H

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
  using namespace TSFExtended;

  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * InserterFactoryBase is a factory interface for element inserters.
     * The create() method builds an inserter of the type associated
     * with this factory. 
     *
     * The create() method is pure virtual and so must be implemented
     * for each derived type. Most inserter factories can be
     * implemented trivially by templating the GenericInserterFactory
     * on the desired inserter type. For example, to make a factory class
     * which creates inserters of type <t>MyInserter</t>, simply 
     * \code
     * RefCountPtr<InserterFactoryBase> factory
     *    = rcp(new GenericInserterFactory<MyInserter>());
     * RefCountPtr<InserterBase> myIns 
     *    = factory.create(rowMap, colMap, bcRows);
     * \endcode
     * 
     */
    class InserterFactoryBase 
    {
    public:
      /** Empty ctor */
      InserterFactoryBase(){;}

      /** Virtual dtor */
      virtual ~InserterFactoryBase(){;}
      
      /** Factory method */
      virtual RefCountPtr<InserterBase> 
      create(const RefCountPtr<DOFMapBase>& rowMap,
             const RefCountPtr<DOFMapBase>& colMap,
             const RefCountPtr<Set<int> >& bcRows) const = 0 ;
                          

    
    };

    /**
     * Generic implementation of the InserterFactoryBase interface.
     */

    template <class T>
    class GenericInserterFactory : public InserterFactoryBase
    {
    public:
      /** */
      GenericInserterFactory() : InserterFactoryBase() {;}

      /** */
      virtual ~GenericInserterFactory(){;}

      /** */
      virtual RefCountPtr<InserterBase> 
      create(const RefCountPtr<DOFMapBase>& rowMap,
             const RefCountPtr<DOFMapBase>& colMap,
             const RefCountPtr<Set<int> >& bcRows) const 
      {
        return rcp(new T(rowMap, colMap, bcRows));
      }
    };


  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
