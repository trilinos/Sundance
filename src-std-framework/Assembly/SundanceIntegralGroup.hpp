/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_INTEGRALGROUP_H
#define SUNDANCE_INTEGRALGROUP_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBuilder.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  namespace Internal
  {
    /** 
     *
     */
    class IntegralGroup
    {
    public:
      /** */
      LocalMatrixContainer(const Array<int>& isTwoForm,
                       const Array<Array<int> >& testID,
                       const Array<Array<int> >& unkID,
                       const Array<Array<double> >& coeffs);

      /** Return the data vector for the i-th batch */
      const RefCountPtr<Array<double> >& dataVector(int i) const 
      {return workspace()[i];}

      /** Indicate whether the i-th batch is a two form */
      bool isTwoForm(int i) const {return isTwoForm_[i];}

      /** Return the array of testIDs whose local matrix values are grouped int
       * the i-th batch */
      const Array<int>& testID(int i) const {return testID_[i];}

      /** Return the array of unkIDs whose local matrix values are grouped int
       * the i-th batch */
      const Array<int>& unkID(int i) const {return unkID_[i];}

      /** Return the array of coefficients to be used with the i-th batch */
      const Array<double>& coeffs(int i) const {return coeffs_[i];}

    private:
      
      /** */
      static Array<RefCountPtr<Array<double> > >& workspace() 
      {static Array<RefCountPtr<Array<double> > > rtn; return rtn;}

      /** */
      Array<int> isTwoForm_;

      /** */
      Array<Array<int> > testID_;

      /** */
      Array<Array<int> > unkID_;

      /** */
      Array<Array<double> > coeffs_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
