/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_INTEGRATORFACTORYBASE_H
#define SUNDANCE_INTEGRATORFACTORYBASE_H

#include "SundanceDefs.hpp"
#include "SundanceIntegratorBase.hpp"

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
     * IntegratorFactoryBase is a factory interface for element integrators.
     * The create() method builds an integrator of the type associated
     * with this factory. 
     *
     * The create() method is pure virtual and so must be implemented
     * for each derived type. Most integrator factories can be
     * implemented trivially by templating the GenericIntegratorFactory
     * on the dersired integrator type. For example, to make a factory class
     * which creates integrators of type <t>MyIntegrator</t>, simply 
     * \code
     * RefCountPtr<IntegratorFactoryBase> factory
     *    = rcp(new GenericIntegratorFactory<MyIntegrator>());
     * RefCountPtr<IntegratorBase> myInt 
     *    = factory.create(mesh,expr,rqc,derivs,mgr);
     * \endcode
     * 
     */
    class IntegratorFactoryBase 
    {
    public:
      /** Empty ctor */
      IntegratorFactoryBase(){;}

      /** Virtual dtor */
      virtual ~IntegratorFactoryBase(){;}
      
      /** Factory method */
      virtual RefCountPtr<IntegratorBase> create(const Mesh& mesh, 
                                                 const Expr& expr,
                                                 const RegionQuadCombo& rqc,
                                                 const DerivSet& derivs,
                                                 const RefCountPtr<EvalManager>& evalMgr) const = 0 ;
                          

    
    };

    /**
     * Generic implementation of the IntegratorFactoryBase interface.
     */

    template <class T>
    class GenericIntegratorFactory : public IntegratorFactoryBase
    {
    public:
      /** */
      GenericIntegratorFactory() : IntegratorFactoryBase() {;}

      /** */
      virtual ~GenericIntegratorFactory(){;}

      /** */
      virtual RefCountPtr<IntegratorBase> create(const Mesh& mesh, 
                                                 const Expr& expr,
                                                 const RegionQuadCombo& rqc,
                                                 const DerivSet& derivs,
                                                 const RefCountPtr<EvalManager>& evalMgr) const 
      {
        return rcp(new T(mesh, expr, derivs, rqc, evalMgr));
      }
    };


  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
