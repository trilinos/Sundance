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

#ifndef SUNDANCE_SPARSITYSUBSET_H
#define SUNDANCE_SPARSITYSUBSET_H


#include "SundanceDefs.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using SundanceUtils::Map;

  namespace Internal
    {
      class SparsitySuperset;
      class EvalVector;

      /**
       * DerivState can be used to classify the known state of
       * each functional derivative at any node in the expression
       * tree.  ZeroDeriv means the derivative is structurally
       * zero at that node -- note that derivatives that are
       * nonzero at a root node can be structurally zero at some
       * child node. ConstantDeriv means that the derivative is
       * known to be a constant (in space) at that
       * node. VectorDeriv means that the derivative is non-zero,
       * non-constant, i.e., a vector of values.
       */
      enum DerivState {ZeroDeriv, ConstantDeriv, VectorDeriv};

      /**
       * Class SparsitySubset contains the subset of functional and spatial
       * derivatives required for evaluation of a particular set of
       * multiindices. 
       */
      class SparsitySubset
        : public TSFExtended::ObjectWithVerbosity<SparsitySubset>
        {
        public:
          /** */
          SparsitySubset(SparsitySuperset* master);

          /** */
          void addDeriv(const MultipleDeriv& d, 
                        const DerivState& state);

          /** */
          void addDeriv(const Deriv& d, 
                        const DerivState& state);
          

          /** Detect whether a given derivative exists in this set */
          bool containsDeriv(const MultipleDeriv& d) const ;

          /** Find the index at which the results for the
              given functional derivative are stored in the results array */
          int getIndex(const MultipleDeriv& d) const ;

          /** Return the derivative stored at index i */
          const MultipleDeriv& deriv(int i) const ;

          /** Return the constancy state of derivative i */
          const DerivState& state(int i) const ;

          /** */
          int numDerivs() const ;

          /** Indicate whether the specified derivative which is
           * nonzero at the root
           * happens to be zero at this node */
          bool isZero(int i) const ;


          /** Indicate whether the specified derivative is
           * spatially constant at this node */
          bool isConstant(int i) const ;

          /** Indicate whether the specified multiple derivative contains
           * at least one order of spatial differentiation */
          bool isSpatialDeriv(int i) const ;

          /** Return the spatial multi index for the i-th derivative */
          const MultiIndex& multiIndex(int i) const ;

          /** */
          void print(ostream& os) const ;

          /** */
          void print(ostream& os, 
                     const Array<RefCountPtr<EvalVector> >& vecResults,
                     const Array<double>& constantResults) const ;

          /** */
          std::string toString() const ;

          /** */
          DerivSet derivSet() const ;

        private:

          const SparsitySuperset* master() const ;

          SparsitySuperset* master() ;

          

          SparsitySuperset* master_;

          Array<MultipleDeriv> derivs_;

          Map<MultipleDeriv, int> derivToIndexMap_;

          Array<int> ptr_;

        };

      


    }
}

namespace std
{
  /** \relates SparsitySubset */
  inline std::ostream& operator<<(std::ostream& os,
                                  const SundanceCore::Internal::SparsitySubset& s)
  {
    s.print(os);
    return os;
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
