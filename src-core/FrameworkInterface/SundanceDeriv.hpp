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

#ifndef SUNDANCE_DERIV_H
#define SUNDANCE_DERIV_H

#include "SundanceDefs.hpp"
#include "SundanceDerivBase.hpp"
#include "SundanceMultiIndex.hpp"
#include "TSFHandle.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY



namespace SundanceCore
{
  using namespace SundanceUtils;

  namespace Internal
    {
      class CoordDeriv;
      class FunctionalDeriv;
      class FuncElementBase;

      using namespace Teuchos;

      /**
       * Deriv represents a first-order functional
       * derivative used in the internal automatic differentiation
       * process.
       * Multiple functional derivatives
       * are represented with the MultipleDeriv object, which is a multiset
       * of Deriv objects.
       *
       * <B>
       * This object should not be confused with the user-level Derivative
       * object, which represents spatial differential operators as appearing
       * in a user-level problem specification.
       * </b>
       *
       * Functional derivatives may be with respect
       * to unknown or test functions and spatial derivatives thereof,
       * or a coordinate function. Under the hood, derivatives with respect
       * to test and unknown functions are represented with
       * FunctionalDeriv objects and derivatives with respect to
       * coordinate functions are represented with CoordDeriv objects.
       * Both FunctionalDeriv and CoordDeriv derive from DerivBase. A Deriv
       * object is a ref-counted handle to a DerivBase.
       *
       * The spatial derivatives represented by Deriv objects are
       * partial derivatives with respect to <b> explicit</b>
       * dependence on a particular coordinate function. For example,
       * the action of a Deriv \f$\mathcal{D}\f$ on
       * a function \f$ F(u,x) \f$ is
       * \f[ \mathcal{D}f = \frac{\partial F}{\partial x}\f]
       * <center><b>not</b></center>
       * \f[ Df = \frac{\partial F}{\partial x}
       * + \frac{\partial F}{\partial u}\frac{\partial u}{\partial x}\f]
       * Chain-rule propagation of implicit dependence on spatial variables
       * is done through automatic differentiation in the evaluation
       * of DiffOp expressions.
       */
      class Deriv : public TSFExtended::Handle<DerivBase>
        {
        public:
          /* handle boilerplate */
          HANDLE_CTORS(Deriv, DerivBase);

          /** */
          bool operator<(const Deriv& other) const ;

          /** */
          bool operator==(const Deriv& other) const ;

          /** */
          string toString() const ;

          /** */
          bool isFunctionalDeriv() const {return funcDeriv() != 0;}

          /** */
          bool isCoordDeriv() const {return coordDeriv() != 0;}

          /** */
          bool isTestFunction() const ;

          /** */
          bool isUnknownFunction() const ;

          /** Return internal pointer downcasted
           * to a FunctionalDeriv pointer */
          const FunctionalDeriv* funcDeriv() const ;

          /** Return internal pointer downcasted to
           * a CoordDeriv pointer */
          const CoordDeriv* coordDeriv() const ;
        };

    }
}

namespace std
{
  /** \relates Internal::Deriv */
  inline ostream& operator<<(ostream& os, 
                             const SundanceCore::Internal::Deriv& d)
    {os << d.toString(); return os;}
}



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
