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
#include "SundanceEvalContext.hpp"
#include "TSFObjectWithVerbosity.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;
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
     * DOF map here because in the Sundance core we know nothing
     * of the mesh, so we provide accessors to the information collected
     * by the EquationSet.
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

      /** \name Methods to access information required for building DOF maps */
      //@{

      /** Returns the number of regions on which pieces of the equation
       * or BCs are defined. */
      int numRegions() const {return regions_.size();}
      
      /** Returns the d-th region for this equation set */
      const RefCountPtr<CellFilterStub>& region(int d) const 
      {return regions_[d].ptr();}

      /** Indicate whether the given region has an essential BC expression */
      bool isBCRegion(int d) const
      {return bcTestsOnRegions_.containsKey(regions_[d]);}
      
      /** Returns the number of test functions in this equation set */
      int numTests() const {return testFuncs_.size();}

      /** Returns the number of unk functions in this equation set */
      int numUnks() const {return unkFuncs_.size();}

      /** Returns the i-th test function */
      const Expr& testFunc(int i) const {return testFuncs_[i];}

      /** Returns the i-th test function */
      const Expr& unkFunc(int i) const {return unkFuncs_[i];}

      /** Returns the test functions that appear explicitly
       * on the d-th region */
      const Set<int>& testsOnRegion(int d) const 
      {return testsOnRegions_.get(regions_[d]);}

      /** Returns the unknown functions that appear explicitly on the
       * d-th region. */
      const Set<int>& unksOnRegion(int d) const 
      {return unksOnRegions_.get(regions_[d]);}

      /** Returns the test functions that appear in BCs on the d-th region.
       * We can use this information to tag certain rows as BC rows */
      const Set<int>& bcTestsOnRegion(int d) const 
      {return bcTestsOnRegions_.get(regions_[d]);}

      /** Returns the unknown functions that appear in BCs on the d-th region.
       * We can use this information to tag certain columns as BC
       * columns in the event we're doing symmetrized BC application */
      const Set<int>& bcUnksOnRegion(int d) const 
      {return bcUnksOnRegions_.get(regions_[d]);}

      /** Determine whether a given func ID is listed as a test function 
       * in this equation set */
      bool hasTestID(int fid) const 
      {return testIDToReducedIDMap_.containsKey(fid);}

      /** Determine whether a given func ID is listed as a unk function 
       * in this equation set */
      bool hasUnkID(int fid) const 
      {return unkIDToReducedIDMap_.containsKey(fid);}

      /** get the reduced ID for the given test ID */
      int reducedTestID(int testID) const 
      {return testIDToReducedIDMap_.get(testID);}

      /** get the reduced ID for the given unk ID */
      int reducedUnkID(int unkID) const 
      {return unkIDToReducedIDMap_.get(unkID);}
      
      //@}


      /** \name Methods needed during setup of integration and evaluation */
      //@{
      /** Returns the list of distinct subregion-quadrature combinations
       * appearing in the equation set. */
      const Array<RegionQuadCombo>& regionQuadCombos() const 
      {return regionQuadCombos_;}

      /** Returns the list of distinct subregion-quadrature combinations
       * appearing in the boundary conditions */
      const Array<RegionQuadCombo>& bcRegionQuadCombos() const 
      {return bcRegionQuadCombos_;}

      /** Returns the set of nonzero functional derivatives appearing
       * in the equation set at the given subregion-quadrature combination */
      const DerivSet& nonzeroFunctionalDerivs(int order,
                                              const RegionQuadCombo& r) const 
      {return regionQuadComboNonzeroDerivs_[order-1].get(r);}

      /** Returns the set of nonzero functional derivatives appearing
       * in the boundary conditions
       *  at the given subregion-quadrature combination */
      const DerivSet& nonzeroBCFunctionalDerivs(int order,
                                                const RegionQuadCombo& r) const
      {return bcRegionQuadComboNonzeroDerivs_[order-1].get(r);}

      /** Map RQC to the context for the derivs of the given order */
      EvalContext rqcToContext(int order, const RegionQuadCombo& r) const 
      {return rqcToContext_[order-1].get(r);}

      /** Map BC RQC to the context for the derivs of the given order */
      EvalContext bcRqcToContext(int order, const RegionQuadCombo& r) const 
      {return bcRqcToContext_[order-1].get(r);}


      /** Indicates whether any test-unk pairs appear in the given domain */
      bool hasTestUnkPairs(const OrderedHandle<CellFilterStub>& domain) const 
      {return testUnkPairsOnRegions_.containsKey(domain);}


      /** Indicates whether any BC test-unk pairs appear in the given domain */
      bool hasBCTestUnkPairs(const OrderedHandle<CellFilterStub>& domain) const 
      {return bcTestUnkPairsOnRegions_.containsKey(domain);}

      /** Returns the (test, unk) pairs appearing on the given domain.
       * This is required for determining the sparsity structure of the
       * matrix */
      const RefCountPtr<Set<OrderedPair<int, int> > >& testUnkPairs(const OrderedHandle<CellFilterStub>& domain) const 
      {return testUnkPairsOnRegions_.get(domain);}
      

       /** Returns the (test, unk) pairs appearing on the given domain.
       * This is required for determining the sparsity structure of the
       * matrix */
      const RefCountPtr<Set<OrderedPair<int, int> > >& bcTestUnkPairs(const OrderedHandle<CellFilterStub>& domain) const ;


      /** */
      const Expr& expr(const RegionQuadCombo& r) const 
      {return regionQuadComboExprs_.get(r);}

      /** */
      const Expr& bcExpr(const RegionQuadCombo& r) const 
      {return bcRegionQuadComboExprs_.get(r);}
      //@}

    private:

      /** */
      void addToTestUnkPairs(const OrderedHandle<CellFilterStub>& domain,
                             const DerivSet& nonzeros, 
                             bool isBC);

      /** */
      Array<OrderedHandle<CellFilterStub> > regions_;

      /** */
      Map<OrderedHandle<CellFilterStub>, Set<int> > testsOnRegions_;

      /** */
      Map<OrderedHandle<CellFilterStub>, Set<int> > unksOnRegions_;

      /** Map from cell filter to pairs of (testID, unkID) appearing
       * on those cells. This is needed to construct the sparsity pattern
       * of the matrix. */
      Map<OrderedHandle<CellFilterStub>, RefCountPtr<Set<OrderedPair<int, int> > > > testUnkPairsOnRegions_;

      /** Map from cell filter to pairs of (testID, unkID) appearing
       * on those cells. This is needed to construct the sparsity pattern
       * of the matrix. */
      Map<OrderedHandle<CellFilterStub>, RefCountPtr<Set<OrderedPair<int, int> > > > bcTestUnkPairsOnRegions_;

      /** */
      Map<OrderedHandle<CellFilterStub>, Set<int> > bcTestsOnRegions_;

      /** */
      Map<OrderedHandle<CellFilterStub>, Set<int> > bcUnksOnRegions_;


      /** */
      Array<RegionQuadCombo> regionQuadCombos_;

      /** */
      Array<RegionQuadCombo> bcRegionQuadCombos_;

      /** */
      Map<RegionQuadCombo, Expr> regionQuadComboExprs_;

      /** */
      Map<RegionQuadCombo, Expr> bcRegionQuadComboExprs_;

      /** List of the sets of nonzero functional derivatives at 
       * each regionQuadCombo */
      Array<Map<RegionQuadCombo, DerivSet> > regionQuadComboNonzeroDerivs_;

      /** List of the sets of nonzero functional derivatives at 
       * each regionQuadCombo */
      Array<Map<RegionQuadCombo, DerivSet> > bcRegionQuadComboNonzeroDerivs_;

      /** List of the contexts for
       * each regionQuadCombo */
      Array<Map<RegionQuadCombo, EvalContext> > rqcToContext_;

      /** List of the contexts for
       * each BC regionQuadCombo */
      Array<Map<RegionQuadCombo, EvalContext> > bcRqcToContext_;

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
