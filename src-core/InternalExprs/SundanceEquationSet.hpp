/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EQUATIONSET_H
#define SUNDANCE_EQUATIONSET_H

#include "SundanceDefs.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSumOfBCs.hpp"
#include "TSFObjectWithVerbosity.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace FrameworkInterface;
  using std::string;

  namespace Internal
    {
      /** 
       *
       */
      class EquationSet : public TSFExtended::ObjectWithVerbosity<EquationSet>
        {
        public:
          /** */
          EquationSet(const Expr& eqns, 
                      const Expr& bcs,
                      const Expr& tests, 
                      const Expr& unks,
                      const Expr& unkLinearizationPts);

          /** */
          const RefCountPtr<CellFilterBase>& domain(int d) const ;

          /** */
          const Array<Expr>& testsOnDomain(int d) const ;

          /** */
          const Array<Expr>& unksOnDomain(int d) const ;

          /** */
          const Array<Expr>& bcTestsOnDomain(int d) const ;

          /** */
          const Array<Expr>& bcUnksOnDomain(int d) const ;
          
        private:
          /** */
          Set<OrderedHandle<CellFilterBase> > domains_;

          /** */
          Map<OrderedHandle<CellFilterBase>, Set<int> > testsOnDomains_;

          /** */
          Map<OrderedHandle<CellFilterBase>, Set<int> > unksOnDomains_;

          /** */
          Map<OrderedHandle<CellFilterBase>, Set<int> > bcTestsOnDomains_;

          /** */
          Map<OrderedHandle<CellFilterBase>, Set<int> > bcUnksOnDomains_;


 //          /** */
//           Array<EvalRegion> regions_;

//           /** */
//           Array<RefCountPtr<CellFilterBase> > regionDomains_;

//           /** */
//           Array<RefCountPtr<QuadratureFamilyBase> > regionQuads_;

//           /** */
//           Array<Exprs> regionExprs_;

//           /** */
//           Array<Exprs> regionBCs_;

//           /** List of the sets of nonzero functional derivatives at 
//            * each region */
//           Array<DerivSet> regionNonzeroDerivs_;

//           /** test functions for this equation set */
//           Expr testFuncs_;

//           /** unknown functions for this equation set */
//           Expr unkFuncs_;

//           /** The point in function space about which the equations
//            * are linearized */
//           Expr unkLinearizationPts_;

//           /** map from test function funcID to that function's
//            * position in list of test functions */
//           Map<int, int> testIDToReducedIDMap_;

//           /** map from unknown function funcID to that function's
//            * position in list of unk functions */
//          Map<int, int> unkIDToReducedIDMap_;
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
