/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_PARAMETER_H
#define SUNDANCE_PARAMETER_H


#include "SundanceDefs.hpp"
#include "SundanceFuncElementBase.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceLeafExpr.hpp"


namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;

  using std::string;
  using std::ostream;

  /** 
   * A Parameter is an expression subtype used to represent
   * a spatially-constant parameter that may change during the
   * course of a simulation, for example, the time in a transient
   * simulation or a continuation parameter when using a homotopy method.
   * While it is possible to use simple double-precision 
   * constants in expressions, their values are immutable once created.
   * When a constant's value may need to be changed, use a Parameter rather
   * than a simple double.
   */
  class Parameter : public Internal::FuncElementBase,
                    public Internal::SpatiallyConstantExpr
  {
  public:
    /** */
    Parameter(const double& value, const string& name="");

    /** virtual destructor */
    virtual ~Parameter() {;}

    /** */
    virtual XMLObject toXML() const ;

#ifndef DOXYGEN_DEVELOPER_ONLY


    /**
     * Indicate whether the given functional derivative is nonzero.
     * A constant expression has a nonzero derivative only if the
     * order of the derivative is zero. 
     */
    virtual bool hasNonzeroDeriv(const MultipleDeriv& d) const
    {return d.order()==0;}

    /**
     * Find all functions and their derivatives beneath my level
     * in the tree. A constant expr has no functions beneath it,
     * so this method does nothing.
     */
    virtual void getRoughDependencies(Set<Deriv>& /* funcs */) const {;}

    /** Write self in text form */
    virtual ostream& toText(ostream& os, bool paren) const 
    {os << "Parameter[" << name() << " = " << value() << "]"; return os;}

    /** */
    virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}

#endif /* DOXYGEN_DEVELOPER_ONLY */
  };
}
#endif
