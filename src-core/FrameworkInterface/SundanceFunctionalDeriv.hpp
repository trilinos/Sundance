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

#ifndef SUNDANCE_FUNCTIONALDERIV_H
#define SUNDANCE_FUNCTIONALDERIV_H

#include "SundanceDefs.hpp"
#include "SundanceDerivBase.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMultiIndex.hpp"
#include "SundanceFuncElementBase.hpp"
#include "Teuchos_RefCountPtr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY



namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      /**
       * FunctionalDeriv represents a single first-order derivative
       * with respect to the derivative of a function, i.e.,
       * \f[\frac{\partial}{\partial \left(D_\alpha u_i\right)}\f]
       * A functional derivative is fully specified with two quantities:
       * the funcID of the function \f$u_i\f$, and the MultiIndex
       * \f$\alpha\f$ which specifies the spatial derivative \f$D_\alpha\f$.
       * These can be accessed by the funcID() and multiIndex() methods,
       * respectively. 
       * @see Deriv
       */
      class FunctionalDeriv : public DerivBase
        {
        public:
          /** Construct with a pointer to a scalar function element and a
           * multiindex */
          FunctionalDeriv(const FuncElementBase* func,
                          const MultiIndex& mi);

          /** virtual dtor */
          virtual ~FunctionalDeriv(){;}

          /** Return the funcID that specifies which function was used in
           * the definition of this functional derivative */
          virtual int funcID() const {return func_->funcID();}

          /** Return the multiindex that specifies which spatial derivative
           *  was used in
           * the definition of this functional derivative*/
          const MultiIndex& multiIndex() const {return mi_;}

          /** write to string */
          virtual string toString() const ;

          /** comparison operator for inclusion in STL ordered containers */
          virtual bool lessThan(const Deriv& other) const ;


          /** Create a new functional derivative in which the argument
           * has been differentiated by the given multi index */
          Deriv derivWrtMultiIndex(const MultiIndex& mi) const ;

          /** */
          const FuncElementBase* func() const {return func_;}

          /* */
          GET_RCP(DerivBase);

        private:
          const FuncElementBase* func_;

          MultiIndex mi_;
        };
    }
}



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
