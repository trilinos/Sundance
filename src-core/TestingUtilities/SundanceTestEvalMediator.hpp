/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TESTEVALMEDIATOR_H
#define SUNDANCE_TESTEVALMEDIATOR_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceTestUnknownFunction.hpp"
#include "SundanceTestDiscreteFunction.hpp"
#include "SundancePoint.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceMap.hpp"

namespace SundanceTesting
{
  using namespace SundanceUtils;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using SundanceUtils::Map;
  /**
   *
   */
  class TestEvalMediator : public AbstractEvalMediator
  {
  public:
    /** */
    TestEvalMediator(const Expr& fields);

    /** */
    virtual ~TestEvalMediator(){;}

    /** */
    void setEvalPoint(const Point& x) {x_=x;}

    /** */
    int numFields() const {return fields_.size();}

    /** */
    void setFieldCoeff(int i, double A) {fields_[i].setCoeff(A);}

    /** */
    double fieldCoeff(int i) const {return fields_[i].coeff();}

    /** */
    const Map<int, int>& funcIdToFieldNumberMap() const
    {return funcIdToFieldNumberMap_;}

    /** Evaluate the given coordinate expression, putting
     * its numerical values in the given LoadableVector. */
    virtual void evalCoordExpr(const CoordExpr* expr,
                               RefCountPtr<EvalVector>& vec) const ;

    /** Evaluate the given cell diameter expression, putting
     * its numerical values in the given EvalVector. */
    virtual void evalCellDiameterExpr(const CellDiameterExpr* expr,
                                      RefCountPtr<EvalVector>& vec) const ;

    /** Evaluate the given discrete function, putting
     * its numerical values in the given LoadableVector. */
    virtual void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                         const Array<MultiIndex>& mi,
                                         Array<RefCountPtr<EvalVector> >& vec) const ;

    double evalDummyBasis(int m, const MultiIndex& mi) const ;


  private:

    Point x_;
    Map<int, int> funcIdToFieldNumberMap_;
    Array<ADField> fields_;
  };
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
