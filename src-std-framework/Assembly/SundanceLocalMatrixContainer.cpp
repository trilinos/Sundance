/* @HEADER@ */
/* @HEADER@ */

#include "SundanceLocalMatrixContainer.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;



LocalMatrixContainer::LocalMatrixContainer(const Array<int>& isTwoForm,
                                           const Array<Array<int> >& testID,
                                           const Array<Array<int> >& unkID,
                                           const Array<Array<double> >& coeffs)
  : isTwoForm_(isTwoForm),
    testID_(testID),
    unkID_(unkID),
    coeffs_(coeffs)
{
  for (int i=workspace().size(); i<isTwoForm.size(); i++)
    {
      workspace().append(rcp(new Array<double>()));
      workspace()[i]->reserve(1000);
    }
}


