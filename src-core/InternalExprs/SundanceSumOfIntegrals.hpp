
/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
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

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;
  using std::string;
  using SundanceUtils::Map;

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

          /** Return the set of unknown or variational
           * functions defined on the d-th region */
          Set<int> funcsOnRegion(int d, const Set<int>& funcsSet) const ;

          /** Indicate whether the integral over the 
           * d-th region contains any test functions */
          bool integralHasTestFunctions(int d) const ;

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
