/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EXPRFIELDWRAPPER_H
#define SUNDANCE_EXPRFIELDWRAPPER_H


#ifndef DOXYGEN_DEVELOPER_ONLY


#include "SundanceDefs.hpp"
#include "TSFHandleable.hpp"
#include "SundanceFieldBase.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceExpr.hpp"

namespace SundanceStdFwk
{
  namespace Internal
  {
    using namespace SundanceCore;
    using namespace SundanceCore::Internal;
    using namespace SundanceStdMesh;
    using namespace SundanceStdMesh::Internal;
    using namespace TSFExtended;
    /**
     *
     */
    class ExprFieldWrapper : public FieldBase
    {
    public:
      /** */
      ExprFieldWrapper(const Expr& expr) ;

      /** virtual dtor */
      virtual ~ExprFieldWrapper(){;}

      /** */
      virtual int numElems() const ;

      /** */
      virtual double getData(int cellDim, int cellID, int elem) const ;

      /* */
      GET_RCP(FieldBase);

    public:
      Expr expr_;

      Vector<double> vector_;

      DiscreteSpace discreteSpace_;

      RefCountPtr<DOFMapBase> map_;
      
      Array<int> indices_;

    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
