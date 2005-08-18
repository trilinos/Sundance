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
     *
     */
    enum ComputationType {MatrixAndVector, VectorOnly, 
                          FunctionalOnly, FunctionalAndGradient};
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
     *
     */
    class EquationSet : public TSFExtended::ObjectWithVerbosity<EquationSet>
    {
    public:
      /** Set up an equation to be integrated with field variables fixed
       * to specified values */
      EquationSet(const Expr& eqns, 
                  const Expr& bcs,
                  const Expr& fields,
                  const Expr& fieldValues);

      /** Set up equations written in weak form with test functions */
      EquationSet(const Expr& eqns, 
                  const Expr& bcs,
                  const Expr& vars, 
                  const Expr& unks,
                  const Expr& unkLinearizationPts);

      /* Set up calculation of a functional and its derivative wrt a 
       * set of variational functions */
      EquationSet(const Expr& eqns, 
                  const Expr& bcs,
                  const Expr& vars,
                  const Expr& varLinearizationPts, 
                  const Expr& fixedFields,
                  const Expr& fixedFieldValues);
      /** Derive a variational problem from a functional */
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
      unsigned int numRegions() const {return regions_.size();}
      
      /** Returns the d-th region for this equation set */
      const RefCountPtr<CellFilterStub>& region(int d) const 
      {return regions_[d].ptr();}

      /** Indicate whether the given region has an essential BC expression */
      bool isBCRegion(int d) const ;
      
      /** Returns the number of variational functions in this equation set */
      unsigned int numVars() const {return varFuncs_.size();}

      /** Returns the number of unk functions in this equation set */
      unsigned int numUnks() const {return unkFuncs_.size();}

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

      /** */
      bool isFunctionalCalculator() const {return isFunctionalCalculator_;}
      

      /** Returns the list of distinct subregion-quadrature combinations
       * appearing in the equation set. */
      const Array<RegionQuadCombo>& regionQuadCombos() const 
      {return regionQuadCombos_;}

      /** Returns the list of distinct subregion-quadrature combinations
       * appearing in the boundary conditions */
      const Array<RegionQuadCombo>& bcRegionQuadCombos() const 
      {return bcRegionQuadCombos_;}

      /** Indicate whether this equation set will do the
       * given computation type */
      bool hasComputationType(ComputationType compType) const 
      {return compTypes_.contains(compType);}

      /** Return the types of computations this object can perform */
      const Set<ComputationType>& computationTypes() const 
      {return compTypes_;}

      /** Returns the set of nonzero functional derivatives appearing
       * in the equation set at the given subregion-quadrature combination */
      const DerivSet& nonzeroFunctionalDerivs(ComputationType compType,
                                              const RegionQuadCombo& r) const ;

      /** Returns the set of nonzero functional derivatives appearing
       * in the boundary conditions
       *  at the given subregion-quadrature combination */
      const DerivSet& nonzeroBCFunctionalDerivs(ComputationType compType,
                                                const RegionQuadCombo& r) const;

      /** Map RQC to the context for the derivs of the given compType */
      EvalContext rqcToContext(ComputationType compType, 
                               const RegionQuadCombo& r) const ;


      /** Map BC RQC to the context for the derivs of the given compType */
      EvalContext bcRqcToContext(ComputationType compType, 
                                 const RegionQuadCombo& r) const ; 


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
                            const Set<int>& vars,
                            const Set<int>& unks,
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
      Map<ComputationType, Map<RegionQuadCombo, DerivSet> > regionQuadComboNonzeroDerivs_;

      /** List of the sets of nonzero functional derivatives at 
       * each regionQuadCombo */
      Map<ComputationType, Map<RegionQuadCombo, DerivSet> > bcRegionQuadComboNonzeroDerivs_;

      /** List of the contexts for
       * each regionQuadCombo */
      Map<ComputationType, Map<RegionQuadCombo, EvalContext> > rqcToContext_;

      /** List of the contexts for
       * each BC regionQuadCombo */
      Map<ComputationType, Map<RegionQuadCombo, EvalContext> > bcRqcToContext_;

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

      /** Set of the computation types supported her */
      Set<ComputationType> compTypes_;

      
      /** Flag indicating whether this equation set is nonlinear */
      bool isNonlinear_;
      
      /** Flag indicating whether this equation set is 
       * a variational problem */
      bool isVariationalProblem_;

      /** Flag indicating whether this equation set is a functional
       * calculator */
      bool isFunctionalCalculator_;

    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
