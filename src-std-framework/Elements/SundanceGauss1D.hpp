/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_GAUSS1D_H
#define SUNDANCE_GAUSS1D_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace Teuchos;
  namespace Internal
  {
    /**
     * Get abscissas and weights for Gauss-Legendre quadrature 
     * on a line segement.
     */
    class Gauss1D
    {
    public:
      /** create an n-point rule on the interval [-1, 1] */
      Gauss1D(int n);
      /** create an n-point rule on the interval [a, b] */
      Gauss1D(int n, double a, double b);

      /** return the number of points in the rule */
      int nPoints() const {return nodes_.length();}
      /** get the abscissas */
      const Array<double>& nodes() const {return nodes_;}
      /** get the weights */
      const Array<double>& weights() const {return weights_;}

      static bool unitTest();
    private:
      void computeWeights(int n, double a, double b);

      Array<double> nodes_;
      Array<double> weights_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
