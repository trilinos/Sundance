/* @HEADER@ */
/* @HEADER@ */


#ifndef SUNDANCE_EVALVECTORARRAY_H
#define SUNDANCE_EVALVECTORARRAY_H


#include "SundanceDefs.hpp"
#include "SundanceEvalVector.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY



namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  namespace Internal
    {
      class SparsitySuperset;

      /**
       *
       */
      class EvalVectorArray : public Array<RefCountPtr<EvalVector> >
        {
        public:
          /** */
          EvalVectorArray();

          /** */
          void copy(const RefCountPtr<EvalVectorArray>& other) ;

          /** */
          void steal(const RefCountPtr<EvalVectorArray>& other);

          /** */
          ostream& print(ostream& os, const SparsitySuperset* derivs) const ;
        };
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif /* EVALRESULTS */

