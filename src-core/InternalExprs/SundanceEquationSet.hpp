/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EQUATIONSET_H
#define SUNDANCE_EQUATIONSET_H

#include "SundanceDefs.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSumOfBCs.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceEvaluatorFactory.hpp"
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
     * EquationSet is an object in which the symbolic specification 
     * of a problem, its BCs, its test and unknown functions, and the
     * point about which it is to be linearized are all gathered. 
     * With this information we can compile lists of which functions
     * are defined on which regions of the domain, which is what is 
     * required for the building of DOF maps. We can't build the
     * DOF map here because in the Sundance core, we know nothing
     * of the mesh, so we must provie accessors to that information.
     */
    class EquationSet : public TSFExtended::ObjectWithVerbosity<EquationSet>
    {
    public:
      /** */
      EquationSet(const Expr& eqns, 
                  const Expr& bcs,
                  const Expr& tests, 
                  const Expr& unks,
                  const Expr& unkLinearizationPts,
                  const RefCountPtr<EvaluatorFactory>& evalFactory);

      /** \name Methods to access information required for building DOF maps */
      //@{

      /** Returns the number of regions on which pieces of the equation
       * or BCs are defined. */
      int numRegions() const ;
      /** Returns the d-th region for this equation set */
      const RefCountPtr<CellFilterBase>& region(int d) const ;

      /** Returns the test functions that appear explicitly
       * on the d-th region */
      const Array<Expr>& testsOnRegion(int d) const ;

      /** Returns the unknown functions that appear explicitly on the
       * d-th region. */
      const Array<Expr>& unksOnRegion(int d) const ;

      /** Returns the test functions that appear in BCs on the d-th region.
       * We can use this information to tag certain rows as BC rows */
      const Array<Expr>& bcTestsOnRegion(int d) const ;

      /** Returns the unknown functions that appear in BCs on the d-th region.
       * We can use this information to tag certain columns as BC
       * columns in the event we're doing symmetrized BC application */
      const Array<Expr>& bcUnksOnRegion(int d) const ;
      //@}

      /** */
      const Array<RegionQuadCombo>& regionQuadCombos() const {return regionQuadComboArray_;}
          
    private:
      /** */
      Set<OrderedHandle<CellFilterBase> > regions_;

      /** */
      Map<OrderedHandle<CellFilterBase>, Set<int> > testsOnRegions_;

      /** */
      Map<OrderedHandle<CellFilterBase>, Set<int> > unksOnRegions_;

      /** */
      Map<OrderedHandle<CellFilterBase>, Set<int> > bcTestsOnRegions_;

      /** */
      Map<OrderedHandle<CellFilterBase>, Set<int> > bcUnksOnRegions_;


      /** */
      Set<RegionQuadCombo> regionQuadCombos_;

      /** */
      Set<RegionQuadCombo> bcRegionQuadCombos_;

      /** */
      Map<RegionQuadCombo, Expr> regionQuadComboExprs_;

      /** */
      Map<RegionQuadCombo, Expr> bcRegionQuadComboExprs_;

      /** List of the sets of nonzero functional derivatives at 
       * each regionQuadCombo */
      Map<RegionQuadCombo, DerivSet> regionQuadComboNonzeroDerivs_;

      /** List of the sets of nonzero functional derivatives at 
       * each regionQuadCombo */
      Map<RegionQuadCombo, DerivSet> bcRegionQuadComboNonzeroDerivs_;

      /** test functions for this equation set */
      Expr testFuncs_;

      /** unknown functions for this equation set */
      Expr unkFuncs_;

      /** The point in function space about which the equations
       * are linearized */
      Expr unkLinearizationPts_;

      /** map from test function funcID to that function's
       * position in list of test functions */
      Map<int, int> testIDToReducedIDMap_;

      /** map from unknown function funcID to that function's
       * position in list of unk functions */
      Map<int, int> unkIDToReducedIDMap_;

      /** Flag indicating whether this equation set is nonlinear */
      bool isNonlinear_;
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
