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

#ifndef SUNDANCE_SPARSITYSUPERSET_H
#define SUNDANCE_SPARSITYSUPERSET_H



#include "SundanceDefs.hpp"
#include "SundanceSparsitySubset.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY



namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  namespace Internal
    {
      

      /**
       * A SparsitySuperset contains a specification of which functional
       * and spatial derivatives are nonzero at a given node in an expression
       * graph in a particular evaluation context. 
       *
       * A superset is further partitioned into subsets containing just 
       * those derivatives required for evaluation of a particular
       * set of spatial derivatives of this expression. 
       */
      class SparsitySuperset 
        : public TSFExtended::ObjectWithVerbosity<SparsitySuperset>
        {
        public:
          
          typedef OrderedPair<Set<MultiIndex>, Set<MultiSet<int> > > keyPair;

          friend class SparsitySubset;

          /** Create an empty superset */
          SparsitySuperset();

          /** \name Subset manipulation */
          //@{
          /** Add a new subset of derivatives, defined as those derivatives
           * required for evaluation of the given set of differential
           * operators applied to this expression */
          void addSubset(const Set<MultiIndex>& multiIndices,
                         const Set<MultiSet<int> >& funcIDs);

          /** Get the subset of derivatives required to evaluate the
           * given set of differential operators */
          const RefCountPtr<SparsitySubset>& subset(const Set<MultiIndex>& multiIndices,
                                                    const Set<MultiSet<int> >& funcIDs) const ;

          /** Get the subset of derivatives required to evaluate the
           * given set of differential operators*/
          RefCountPtr<SparsitySubset> subset(const Set<MultiIndex>& multiIndices,
                                             const Set<MultiSet<int> >& funcIDs) ;

          /** Tell whether the specified subset has been defined */
          bool hasSubset(const Set<MultiIndex>& multiIndices,
                         const Set<MultiSet<int> >& funcIDs) const ;
          //@}

          /** \name Access to information about individual derivatives */
          //@{
          /** Detect whether a given derivative exists in this set */
          bool containsDeriv(const MultipleDeriv& d) const ;

          /** Find the index at which the results for the
              given functional derivative are stored in the results array */
          int getIndex(const MultipleDeriv& d) const ;

          /** Return the results stored at index i */
          inline const MultipleDeriv& deriv(int i) const
            {return derivs_[i];}

          /** Return the constancy state of deriv i */
          inline const DerivState& state(int i) const
            {return states_[i];}

          /** */
          inline int numDerivs() const {return derivs_.size();}

          /** */
          int numConstantDerivs() const {return numConstantDerivs_;}

          /** */
          int numVectorDerivs() const {return numVectorDerivs_;}

          
          
          /** */
          int maxOrder() const {return maxOrder_;}


          /** Indicate whether the specified derivative is
           * spatially constant at this node */
          bool isConstant(int i) const {return states_[i]==ConstantDeriv;}

          /** Indicate whether the specified multiple derivative contains
           * at least one order of spatial differentiation */
          bool isSpatialDeriv(int i) const 
          {
            return multiIndex_[i].order() != 0;
          }

          /** Return the spatial multi index for the i-th derivative */
          const MultiIndex& multiIndex(int i) const 
          {return multiIndex_[i];}
          //@}
          
          const Set<MultiIndex>& allMultiIndices() const 
          {return allMultiIndices_;}

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

          /** */
          void displayAll(ostream& os) const ;

        private:

          /** Add a derivative to the set. Called by the subset's addDeriv()
           * method. */
          void addDeriv(const MultipleDeriv& d, 
                        const DerivState& state);

          /** Add a derivative to the set. Called by the subset's addDeriv()
           * method. */
          void addDeriv(const Deriv& d, 
                        const DerivState& state);

          /** */
          int maxOrder_;

          /** Map from deriv to position of the derivative's
           * value in the results array */
          Map<MultipleDeriv, int> derivToIndexMap_;

          /** The list of functional derivatives whose values are
           * stored in this results set */
          Array<MultipleDeriv> derivs_;

          /** The state of each derivative at this node in the expression */
          Array<DerivState> states_;

          /** Multiindices */
          Array<MultiIndex> multiIndex_;

          /** Table of subsets */
          Map<keyPair, RefCountPtr<SparsitySubset> > subsets_;

          /** */
          Set<MultiIndex> allMultiIndices_;

          /** */
          Set<MultiSet<int> > allFuncIDs_;

          /** */
          int numConstantDerivs_;

          int numVectorDerivs_;

          

        };
  }
}

namespace std
{
  /** \relates SparsitySuperset */
  inline std::ostream& operator<<(std::ostream& os,
                                  const SundanceCore::Internal::SparsitySuperset& s)
  {
    s.print(os);
    return os;
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
