/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_INTEGRALGROUP_H
#define SUNDANCE_INTEGRALGROUP_H

#include "SundanceDefs.hpp"
#include "SundanceElementIntegral.hpp"
#include "SundanceEvalVectorArray.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  namespace Internal
  {
    /** 
     *
     */
    class IntegralGroup : public TSFExtended::ObjectWithVerbosity<IntegralGroup>
    {
    public:
      /** */
      IntegralGroup(const Array<int>& testID,
                    const Array<RefCountPtr<ElementIntegral> >& integrals,
                    const Array<Array<int> >& resultIndices);
      /** */
      IntegralGroup(const Array<int>& testID,
                    const Array<int>& unkID,
                    const Array<RefCountPtr<ElementIntegral> >& integrals,
                    const Array<Array<int> >& resultIndices);


      /** Indicate whether this is a group of two-forms */
      bool isTwoForm() const {return isTwoForm_;}

      /** Return the number of rows in the local matrices or vectors
       * computed by this integral group */
      int nTestNodes() const {return nTestNodes_;}

      /** Return the number of columns in the local matrices 
       * computed by this integral group */
      int nUnkNodes() const {return nUnkNodes_;}

      /** Return the test functions using this integral group */
      const Array<int>& testID() const {return testID_;}

      /** Return the unknown functions using this integral group */
      const Array<int>& unkID() const {return unkID_;}

      /** Evaluate this integral group */
      bool evaluate(const CellJacobianBatch& J,
                    const RefCountPtr<EvalVectorArray>& coeffs,
                    RefCountPtr<Array<double> >& A) const ;


    private:
      
      /** */
      bool isTwoForm_;

      /** */
      int nTestNodes_;

      /** */
      int nUnkNodes_;

      /** */
      Array<int> testID_;

      /** */
      Array<int> unkID_;

      /** */
      Array<RefCountPtr<ElementIntegral> > integrals_;

      /** */
      Array<Array<int> > resultIndices_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif