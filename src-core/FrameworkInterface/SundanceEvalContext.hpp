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

#ifndef SUNDANCE_EVALCONTEXT_H
#define SUNDANCE_EVALCONTEXT_H


#include "SundanceDefs.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_Utils.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY



namespace SundanceCore
{
  using namespace Teuchos;
  using namespace SundanceUtils;
  using std::string;

  namespace Internal
    {
      /** 
       * Different contexts might require the same expression to be
       * evaluated to different orders of functional differentiation; for
       * example, in setting up a linear system, second-order derivatives
       * are required, but in evaluating a functional only zeroth derivs
       * are required. 
       * An EvaluationContext is used as a key to associate an evaluator and
       * its corresponding set of
       * functional derivatives with a context.
       *
       * They key consists of three parts: first, an integer identifier
       * indicating the caller, e.g., an assembler or functional evaluator,
       * second, an integer representing the maximum order of 
       * differentiation required by the top level caller, and third,
       a region-quadrature combination.  
       */
      class EvalContext
        {
        public:
          /** Empty ctor */
          EvalContext() : data_() {;}

          /** Construct with a region-quadrature combination and
           * an identifier of the construcing context. */
          EvalContext(const RegionQuadCombo& rqc,
                      int topLevelDiffOrder,
                      int contextID)
            : data_(rcp(new OrderedTriple<int, int, RegionQuadCombo>(topLevelDiffOrder, contextID, rqc)))
          {;}

          /** Comparison operator for use in maps */
          bool operator<(const EvalContext& other) const 
          {return *data_ < *other.data_;}
          
          /** Write to a string */
          string toString() const
          {return "EvalContext[diffOrder=" 
             + Teuchos::toString(data_->a())
             + ", id=" 
             + Teuchos::toString(data_->b())
             + ", " + data_->c().toString() + "]";}
          
          /** Write a short description to a string */
          string brief() const
          {return "EvalContext[diffOrder=" 
             + Teuchos::toString(data_->a())
             + ", id=" 
             + Teuchos::toString(data_->b())
             + "]";}

          /** */
          int topLevelDiffOrder() const {return data_->a();}

          /** Return a unique context ID */
          static int nextID() {static int rtn=0; return rtn++;}
        private:
          RefCountPtr<OrderedTriple<int, int, RegionQuadCombo> > data_;
        };

    }
}


namespace std
{
  /** \relates SundanceCore::Internal::EvalContext */
  inline ostream& operator<<(ostream& os, 
                             const SundanceCore::Internal::EvalContext& c)
  {
    os << c.toString();
    return os;
  }
}

namespace Teuchos
{
  using std::string;

  /** \relates SundanceCore::Internal::EvalContext */
  inline string toString(const SundanceCore::Internal::EvalContext& h)
    {return h.toString();}

}



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
