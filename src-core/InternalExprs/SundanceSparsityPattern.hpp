/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SPARSITYPATTERN_H
#define SUNDANCE_SPARSITYPATTERN_H



#include "SundanceDefs.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY





namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  namespace Internal
    {
      class EvaluatableExpr;

      /**
       * Class SparsityPattern specifies which functional derivatives
       * exist at each node.
       */
      class SparsityPattern
        {
        public:
          /** Given a set of derivatives that should be evaluated at
           * the given expr, construct a sparsity pattern */
          //          SparsityPattern(const DerivSet& d,
          //                          const EvaluatableExpr* expr);

          SparsityPattern(const DerivSet& d, const EvaluatableExpr* e);

          /** Detect whether a given derivative exists in this set */
          bool containsDeriv(const MultipleDeriv& d) const ;

          /** Find the index at which the results for the
              given functional derivative are stored in the results array */
          int getIndex(const MultipleDeriv& d) const ;

          /** Return the results stored at index i */
          inline const MultipleDeriv& deriv(int i) const
            {return derivs_[i];}

          inline int numDerivs() const {return derivs_.size();}

          /** Indicate whether the specified derivative which is
           * nonzero at the root
           * happens to be zero at this node */
          bool isZero(int i) const {return states_[i]==ZeroDeriv;}


          /** Indicate whether the specified derivative is
           * spatially constant at this node */
          bool isConstant(int i) const {return states_[i]==ConstantDeriv;}

          /**
           * Indicate whether the specified derivative is 
           * a first-order spatial deriv.  
           */
          bool isFirstOrderSpatialDeriv(int i) const 
          {
            return isFirstOrderSpatialDeriv_[i];
          }

          /**
           * Return the direction of spatial differentiation.
           */
          int spatialDerivDir(int i) const 
          {
            return spatialDerivDir_[i];
          }


          /** */
          void print(ostream& os) const ;

          /** */
          std::string toString() const ;

        private:

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


          /** Map from deriv to position of the derivative's
           * value in the results array */
          Map<MultipleDeriv, int> derivToIndexMap_;

          /** The list of functional derivatives whose values are
           * stored in this results set */
          Array<MultipleDeriv> derivs_;

          /** The state of each derivative at this node in the expression */
          Array<DerivState> states_;

          /** Flags indicating whether a given derivative is first order wrt 
           * a spatial coordinate */
          Array<int> isFirstOrderSpatialDeriv_;

          /** Directions of spatial derivatives */
          Array<int> spatialDerivDir_;

        };

      /** \relates SparsityPattern */
      inline std::ostream& operator<<(std::ostream& os,
                                      const SparsityPattern& s)
        {
          s.print(os);
          return os;
        }


    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
