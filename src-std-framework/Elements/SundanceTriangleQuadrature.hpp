/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TRIANGLEQUADRATURE_H
#define SUNDANCE_TRIANGLEQUADRATURE_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace Teuchos;
  namespace Internal
  {
    /**
     * Get abscissas and weights for Gaussian quadrature on triangles
     */

    class TriangleQuadrature
    {
    public:
      static void getPoints(int order, Array<double>& wgt,
                            Array<double>& x,
                            Array<double>& y);

      static bool test(int p);


    private:

      static void permute(int m, const Array<double>& q,
                          Array<Array<double> >& qPerm);

      static double exact(int a, int b, int c);

      static double fact(int x);

    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
