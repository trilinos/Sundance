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
     *
     * There are several modes in which one might construct an equation set.
     * The first is where one has written out a weak form in terms
     * of test functions. The second is where one is taking variations
     * of some functional.
     */
    class EquationSet : public TSFExtended::ObjectWithVerbosity<EquationSet>
    {
    public:
      /** */
      EquationSet(const Expr& eqns, 
                  const Expr& bcs,
                  const Expr& vars, 
                  const Expr& unks,
                  const Expr& unkLinearizationPts);
      /** */
      EquationSet(const Expr& eqns, 
                  const Expr& bcs,
                  const Expr& vars, 
                  const Expr& varLinearizationPts,
                  const Expr& unks,
                  const Expr& unkLinearizationPts,
                  const Expr& fixedFields,
                  const Expr& fixedFieldValues);

      /** \name Methods to access information required for building DOF maps */
      //@{

      /** Returns the number of regions on which pieces of the equation
       * or BCs are defined. */
      int numRegions() const {return regions_.size();}
      
      /** Returns the d-th region for this equation set */
      const RefCountPtr<CellFilterStub>& region(int d) const 
      {return regions_[d].ptr();}

      /** Indicate whether the given region has an essential BC expression */
      bool isBCRegion(int d) const ;
      
      /** Returns the number of variational functions in this equation set */
      int numVars() const {return varFuncs_.size();}

      /** Returns the number of unk functions in this equation set */
      int numUnks() const {return unkFuncs_.size();}

      /** Returns the i-th variational function */
      const Expr& varFunc(int i) const {return varFuncs_[i];}

      /** Returns the i-th unknown function */
      const Expr& unkFunc(int i) const {return unkFuncs_[i];}

      /** Returns the variational functions that appear explicitly
       * on the d-th region */
      const Set<int>& varsOnRegion(int d) const 
      {return varsOnRegions_.get(regions_[d]);}

      /** Returns the unknown functions that appear explicitly on the
       * d-th region. */
      const Set<int>& unksOnRegion(int d) const 
      {return unksOnRegions_.get(regions_[d]);}

      /** Returns the variational functions that 
       * appear in BCs on the d-th region.
       * We can use this information to tag certain rows as BC rows */
      const Set<int>& bcVarsOnRegion(int d) const 
      {return bcVarsOnRegions_.get(regions_[d]);}

      /** Returns the unknown functions that appear in BCs on the d-th region.
       * We can use this information to tag certain columns as BC
       * columns in the event we're doing symmetrized BC application */
      const Set<int>& bcUnksOnRegion(int d) const 
      {return bcUnksOnRegions_.get(regions_[d]);}

      /** Determine whether a given func ID is listed as a 
       * variational function in this equation set */
      bool hasVarID(int fid) const 
      {return varIDToReducedIDMap_.containsKey(fid);}

      /** Determine whether a given func ID is listed as a unk function 
       * in this equation set */
      bool hasUnkID(int fid) const 
      {return unkIDToReducedIDMap_.containsKey(fid);}

      /** get the reduced ID for the given variational ID */
      int reducedVarID(int varID) const 
      {return varIDToReducedIDMap_.get(varID);}

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


      /** Indicates whether any var-unk pairs appear in the given domain */
      bool hasVarUnkPairs(const OrderedHandle<CellFilterStub>& domain) const 
      {return varUnkPairsOnRegions_.containsKey(domain);}


      /** Indicates whether any BC var-unk pairs appear in the given domain */
      bool hasBCVarUnkPairs(const OrderedHandle<CellFilterStub>& domain) const 
      {return bcVarUnkPairsOnRegions_.containsKey(domain);}

      /** Returns the (var, unk) pairs appearing on the given domain.
       * This is required for determining the sparsity structure of the
       * matrix */
      const RefCountPtr<Set<OrderedPair<int, int> > >& varUnkPairs(const OrderedHandle<CellFilterStub>& domain) const 
      {return varUnkPairsOnRegions_.get(domain);}
      

       /** Returns the (var, unk) pairs appearing on the given domain.
       * This is required for determining the sparsity structure of the
       * matrix */
      const RefCountPtr<Set<OrderedPair<int, int> > >& bcVarUnkPairs(const OrderedHandle<CellFilterStub>& domain) const ;


      /** */
      const Expr& expr(const RegionQuadCombo& r) const 
      {return regionQuadComboExprs_.get(r);}

      /** */
      const Expr& bcExpr(const RegionQuadCombo& r) const 
      {return bcRegionQuadComboExprs_.get(r);}
      //@}

    private:

      /** */
      void init(const Expr& eqns, 
                const Expr& bcs,
                const Expr& vars, 
                const Expr& varLinearizationPts,
                const Expr& unks,
                const Expr& unkLinearizationPts,
                const Expr& fixedFields,
                const Expr& fixedFieldValues);

      /** */
      void addToVarUnkPairs(const OrderedHandle<CellFilterStub>& domain,
                             const DerivSet& nonzeros, 
                             bool isBC);

      /** */
      Array<OrderedHandle<CellFilterStub> > regions_;

      /** */
      Map<OrderedHandle<CellFilterStub>, Set<int> > varsOnRegions_;

      /** */
      Map<OrderedHandle<CellFilterStub>, Set<int> > unksOnRegions_;

      /** Map from cell filter to pairs of (varID, unkID) appearing
       * on those cells. This is needed to construct the sparsity pattern
       * of the matrix. */
      Map<OrderedHandle<CellFilterStub>, RefCountPtr<Set<OrderedPair<int, int> > > > varUnkPairsOnRegions_;

      /** Map from cell filter to pairs of (varID, unkID) appearing
       * on those cells. This is needed to construct the sparsity pattern
       * of the matrix. */
      Map<OrderedHandle<CellFilterStub>, RefCountPtr<Set<OrderedPair<int, int> > > > bcVarUnkPairsOnRegions_;

      /** */
      Map<OrderedHandle<CellFilterStub>, Set<int> > bcVarsOnRegions_;

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

      /** var functions for this equation set */
      Expr varFuncs_;

      /** unknown functions for this equation set */
      Expr unkFuncs_;

      /** The point in function space about which the equations
       * are linearized */
      Expr unkLinearizationPts_;

      /** map from variational function funcID to that function's
       * position in list of var functions */
      Map<int, int> varIDToReducedIDMap_;

      /** map from unknown function funcID to that function's
       * position in list of unk functions */
      Map<int, int> unkIDToReducedIDMap_;

      
      /** Flag indicating whether this equation set is nonlinear */
      bool isNonlinear_;
      
      /** Flag indicating whether this equation set is 
       * a variational problem */
      bool isVariationalProblem_;
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
