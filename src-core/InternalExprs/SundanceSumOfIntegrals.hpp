/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SUMOFINTEGRALS_H
#define SUNDANCE_SUMOFINTEGRALS_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundanceOrderedHandle.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceMap.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY
using SundanceUtils::Map;
namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;
  using std::string;

  namespace Internal
    {
      class SpatiallyConstantExpr;

      /** 
       * SumOfIntegrals represents a sum of integrals,
       * grouped by region and quadrature rule
       *
       * \f[
       * \sum_{d=0}^{N_d-1} \left[\sum_{q=0}^{N_{q,d}-1} 
       * \int_{\Omega_d,Q_{q,d}} g_{d,q}\right] 
       * \f] 
       *
       * 
       * 
       */
      class SumOfIntegrals : public ScalarExpr
        {
        public:
          /** Construct given an integral over a single region */
          SumOfIntegrals(const RefCountPtr<CellFilterStub>& region,
                         const Expr& expr,
                         const RefCountPtr<QuadratureFamilyStub>& quad);

          /** */
          virtual ~SumOfIntegrals(){;}

          /** Add another term to this integral */
          void addTerm(const RefCountPtr<CellFilterStub>& region,
                       const Expr& expr,
                       const RefCountPtr<QuadratureFamilyStub>& quad, 
                       int sign) ;

          /** Add this sum of integrals to another sum of integrals */
          void merge(const SumOfIntegrals* other, int sign) ;

          /** Multiply all terms in the sum by a constant */
          void multiplyByConstant(const SpatiallyConstantExpr* expr) ;

          /** Change the sign of all terms in the sum */
          void changeSign() ;

          /** Return the number of subregions */
          int numRegions() const {return regions_.size();}

          /** Return the d-th region */
          const RefCountPtr<CellFilterStub>& region(int d) const 
          {return regions_[d].ptr();}

          /** Return the number of different quadrature rules used in the
           * d-th region */
          int numTerms(int d) const {return quad_[d].size();}

          /** Return the q-th quadrature rule used in the d-th region */
          const RefCountPtr<QuadratureFamilyStub>& quad(int d, int q) const 
          {return quad_[d][q].ptr();}

          /** Return the integrand for the q-th quadrature rule in the
           * d-th region */
          const Expr& expr(int d, int q) const 
          {return expr_[d][q];}

          /** Return the set of unknown functions defined on the d-th region */
          Set<int> unksOnRegion(int d) const ;

          /** Return the set of test functions defined on the d-th region */
          Set<int> testsOnRegion(int d) const ;

          /** Return a null cell filter of a type consistent with the
           * other filters in this integral */
          RefCountPtr<CellFilterStub> nullRegion() const ;

          /** Write a simple text description suitable 
           * for output to a terminal */
          virtual ostream& toText(ostream& os, bool paren) const ;

          /** Write in a form suitable for LaTeX formatting */
          virtual ostream& toLatex(ostream& os, bool paren) const ;

          /** Write in XML */
          virtual XMLObject toXML() const ;

          /** */
          virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

        private:
          Array<OrderedHandle<CellFilterStub> > regions_;

          Array<Array<OrderedHandle<QuadratureFamilyStub> > > quad_;

          Array<Array<Expr> > expr_;

          Map<OrderedHandle<CellFilterStub>, int>  cellSetToIndexMap_;

          Array<Map<OrderedHandle<QuadratureFamilyStub>, int> > quadToIndexMap_;
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
