/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALCONTEXT_H
#define SUNDANCE_EVALCONTEXT_H


#include "SundanceDefs.hpp"
#include "SundanceRegionQuadCombo.hpp"
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
            : data_(rcp(new OrderedTriple<int, int, RegionQuadCombo>(topLevelDiffOrder, contextID, rqc))) {;}

          /** Comparison operator for use in maps */
          bool operator<(const EvalContext& other) const 
          {return *data_ < *other.data_;}
          
          /** Write to a string */
          string toString() const
          {return "EvalContext[diffOrder=" 
             + Teuchos::toString(data_->a())
             + "id=" 
             + Teuchos::toString(data_->b())
             + ", " + data_->c().toString();}

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
