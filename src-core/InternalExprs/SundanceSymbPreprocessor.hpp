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



#ifndef SUNDANCE_SYMBPREPROCESSOR_H
#define SUNDANCE_SYMBPREPROCESSOR_H

#include "SundanceExpr.hpp"
#include "SundanceEvaluatableExpr.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;

  using std::string;
  using std::ostream;

  namespace Internal
    {
      /** */
      class SymbPreprocessor
      {
      public:

        /** */
        static DerivSet setupVariations(const Expr& expr, 
                                        const Expr& vars,
                                        const Expr& varEvalPts,
                                        const Expr& unks,
                                        const Expr& unkEvalPts,
                                        const Expr& fixedFields,
                                        const Expr& fixedFieldEvalPts,
                                        const EvalContext& region);

        /** */
        static DerivSet setupFunctional(const Expr& expr, 
                                        const Expr& fixedFields,
                                        const Expr& fixedFieldEvalPts,
                                        const EvalContext& region);

        /** */
        static DerivSet setupGradient(const Expr& expr, 
                                      const Expr& vars,
                                      const Expr& varEvalPts,
                                      const Expr& fixedFields,
                                      const Expr& fixedFieldEvalPts,
                                      const EvalContext& region);

        /** */
        static DerivSet setupExpr(const Expr& expr, 
                                  const Expr& tests,
                                  const Expr& unks,
                                  const Expr& u0,
                                  const EvalContext& region);
        /** */
        static DerivSet setupExpr(const Expr& expr, 
                                  const Expr& unks,
                                  const Expr& u0,
                                  const EvalContext& region);
        /** */
        static DerivSet setupExpr(const Expr& expr, 
                                  const EvalContext& region);

        /** */
        TEUCHOS_TIMER(preprocTimer, "symbolic preprocessing");
      };
  }
}



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
