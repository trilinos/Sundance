/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_COORDEXPR_H
#define SUNDANCE_COORDEXPR_H

#include "SundanceFuncElementBase.hpp"
#include "SundanceLeafExpr.hpp"

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace SundanceCore::Internal;

  /** */
  class CoordExpr : public FuncElementBase,
                    virtual public LeafExpr
    {
    public:
      /** */
      CoordExpr(int dir, const string& name="");

      /** */
      virtual ~CoordExpr() {;}

      /** */
      virtual XMLObject toXML() const ;

#ifndef DOXYGEN_DEVELOPER_ONLY
      /** */
      int dir() const {return dir_;}


      /**
       * Indicate whether the given functional derivative is nonzero
       */
      virtual bool hasNonzeroDeriv(const MultipleDeriv& f) const ;

      /**
       * Find all functions and their derivatives beneath my level
       * in the tree.
       */
      virtual void getRoughDependencies(Set<Deriv>& funcs) const {;}



      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

    private:
      int dir_;

      static string coordName(int dir, const string& name);
#endif  /* DOXYGEN_DEVELOPER_ONLY */

    };
}

#endif
