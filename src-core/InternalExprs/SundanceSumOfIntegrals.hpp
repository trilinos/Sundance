/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SUMOFINTEGRALS_H
#define SUNDANCE_SUMOFINTEGRALS_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "SundanceCellFilterBase.hpp"
#include "SundanceOrderedHandle.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceMap.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace FrameworkInterface;
  using std::string;

  namespace Internal
    {
      class SpatiallyConstantExpr;

      /** 
       * SumOfIntegrals represents a sum of integrals,
       * grouped by domain and quadrature rule
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
          /** Construct given an integral over a single domain */
          SumOfIntegrals(const RefCountPtr<CellFilterBase>& domain,
                         const Expr& expr,
                         const RefCountPtr<QuadratureFamilyBase>& quad);

          /** */
          virtual ~SumOfIntegrals(){;}

          /** Add another term to this integral */
          void addTerm(const RefCountPtr<CellFilterBase>& domain,
                       const Expr& expr,
                       const RefCountPtr<QuadratureFamilyBase>& quad, 
                       int sign) ;

          /** Add this sum of integrals to another sum of integrals */
          void merge(const SumOfIntegrals* other, int sign) ;

          /** Multiply all terms in the sum by a constant */
          void multiplyByConstant(const SpatiallyConstantExpr* expr) ;

          /** Change the sign of all terms in the sum */
          void changeSign() ;

          /** Return the number of subdomains */
          int numDomains() const {return domains_.size();}

          /** Return the d-th domain */
          const RefCountPtr<CellFilterBase>& domain(int d) const 
          {return domains_[d].ptr();}

          /** Return the number of different quadrature rules used in the
           * d-th domain */
          int numTerms(int d) const {return quad_[d].size();}

          /** Return the q-th quadrature rule used in the d-th domain */
          const RefCountPtr<QuadratureFamilyBase>& quad(int d, int q) const 
          {return quad_[d][q].ptr();}

          /** Return the integrand for the q-th quadrature rule in the
           * d-th domain */
          const Expr& expr(int d, int q) const 
          {return expr_[d][q];}

          /** Return the set of unknown functions defined on the d-th domain */
          Set<int> unksOnDomain(int d) const ;

          /** Return the set of test functions defined on the d-th domain */
          Set<int> testsOnDomain(int d) const ;

          /** Return a null cell filter of a type consistent with the
           * other filters in this integral */
          RefCountPtr<CellFilterBase> nullDomain() const ;

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
          Array<OrderedHandle<CellFilterBase> > domains_;

          Array<Array<OrderedHandle<QuadratureFamilyBase> > > quad_;

          Array<Array<Expr> > expr_;

          Map<OrderedHandle<CellFilterBase>, int>  cellSetToIndexMap_;

          Array<Map<OrderedHandle<QuadratureFamilyBase>, int> > quadToIndexMap_;
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
