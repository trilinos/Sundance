/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_WEAKFORMBATCH_H
#define SUNDANCE_WEAKFORMBATCH_H

#include "SundanceDefs.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceBasisFamily.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "TSFPrintable.hpp"

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
    class WeakFormBatch : public TSFExtended::ObjectWithVerbosity<WeakFormBatch>,
                          public TSFExtended::Printable
    {
    public:
      /** */
      WeakFormBatch(const MultipleDeriv& functionalDeriv, int derivIndex);

      /** */
      bool tryToAdd(const MultipleDeriv& functionalDeriv, int derivIndex);

      /** */
      bool isOneForm() const {return isOneForm_;}

      /** */
      const BasisFamily& testBasis() const {return testBasis_;}

      /** */
      const MultiIndex& miTest() const {return miTest_;}

      /** */
      const Array<int>& testID() const {return testID_;}

      /** */
      const BasisFamily& unkBasis() const {return unkBasis_;}

      /** */
      const MultiIndex& miUnk() const {return miUnk_;}

      /** */
      const Array<int>& unkID() const {return unkID_;}

      /** */
      const Array<int>& derivIndex() const {return derivIndex_;}

      /** */
      void print(ostream& os) const ;

    private:
      void extractWeakForm(const MultipleDeriv& functionalDeriv,
                           BasisFamily& testBasis, BasisFamily& unkBasis,
                           MultiIndex& miTest, MultiIndex& miUnk,
                           int& testID, int& unkID, bool& isOneForm) const ;
      bool isOneForm_;
      
      BasisFamily testBasis_;

      BasisFamily unkBasis_;

      MultiIndex miTest_;

      MultiIndex miUnk_;

      Array<int> testID_;

      Array<int> unkID_;

      Array<int> derivIndex_;

      Array<int> trans_;
      
    };

    /** */
    inline ostream& operator<<(ostream& os, const WeakFormBatch& w)
    {
      w.print(os);
      return os;
    }
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
