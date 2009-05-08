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

#ifndef SUNDANCE_ASSEMBLER_H
#define SUNDANCE_ASSEMBLER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceGrouperBase.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "TSFLoadableVector.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFVector.hpp"
#include "TSFVectorType.hpp"
#include "Teuchos_HashSet.hpp"
#include "Teuchos_ParameterList.hpp"
#include "TSFIncrementallyConfigurableMatrixFactory.hpp"
#include "TSFCollectivelyConfigurableMatrixFactory.hpp"
#include "TSFPartitionedMatrixFactory.hpp"
#include "TSFPartitionedToMonolithicConverter.hpp"
#include "SundanceAssemblyKernelBase.hpp"

namespace SundanceStdFwk
{
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceCore;
using namespace SundanceCore;

namespace Internal
{
using namespace Teuchos;

typedef std::set<int> ColSetType;

/** 
 * 
 */
class Assembler 
  : public TSFExtended::ParameterControlledObjectWithVerbosity<Assembler>
{
public:
  /** */
  Assembler(
    const Mesh& mesh, 
    const RefCountPtr<EquationSet>& eqn,
    const Array<VectorType<double> >& rowVectorType,
    const Array<VectorType<double> >& colVectorType,
    bool partitionBCs,
    const ParameterList& verbParams = *defaultVerbParams());


  /** */
  Assembler(
    const Mesh& mesh, 
    const RefCountPtr<EquationSet>& eqn,
    const ParameterList& verbParams= *defaultVerbParams());
      
  /** */
  const Array<RefCountPtr<DOFMapBase> >& rowMap() const 
    {return rowMap_;}

  /** */
  const Array<RefCountPtr<DOFMapBase> >& colMap() const 
    {return colMap_;}

  /** */
  const Array<RefCountPtr<DiscreteSpace> >& solutionSpace() const 
    {return externalColSpace_;}

  /** */
  const Array<RefCountPtr<DiscreteSpace> >& rowSpace() const 
    {return externalRowSpace_;}

  /** */
  VectorSpace<double> solnVecSpace() const ;

  /** */
  VectorSpace<double> rowVecSpace() const ;

  /** */
  const Array<RefCountPtr<Set<int> > >& bcRows() {return bcRows_;}

  /** Allocate, but do not fill, the matrix */
  TSFExtended::LinearOperator<double> allocateMatrix() const ;

  /** */
  void assemble(TSFExtended::LinearOperator<double>& A,
    TSFExtended::Vector<double>& b) const ;


  /** */
  void assemble(TSFExtended::Vector<double>& b) const ;

  /** */
  void evaluate(double& value,
    TSFExtended::Vector<double>& gradient) const ;

  /** */
  void evaluate(double& value) const ;

  /** */
  static unsigned int& workSetSize() ;

      
  /** */
  void getGraph(int br, int bc,
    Array<int>& graphData,
    Array<int>& rowPtrs,
    Array<int>& nnzPerRow) const ;
      
  /** */
  void incrementalGetGraph(int br, int bc, 
    IncrementallyConfigurableMatrixFactory* mf) const ;

  /** */
  void flushConfiguration() 
    {
      vecNeedsConfiguration_ = true;
      matNeedsConfiguration_ = true;
    }

  /** */
  Vector<double> convertToMonolithicVector(
    const Array<Vector<double> >& internalBlock,
    const Array<Vector<double> >& bcBlock) const ;

  /** */
  static int& numAssembleCalls() {static int rtn=0; return rtn;}

  /** */
  static bool& matrixEliminatesRepeatedCols() {static bool x = false; return x;}

  /** */
  const RefCountPtr<EquationSet>& eqnSet() const 
    {return eqn_;}

  /** */
  static RefCountPtr<ParameterList> defaultVerbParams()
    {
      static RefCountPtr<ParameterList> rtn = rcp(new ParameterList("Assembler"));
      static int first = true;
      if (first)
      {
        rtn->set<int>("setup", 0);
        rtn->set<int>("matrix config", 0);
        rtn->set<int>("vector config", 0);
        rtn->set<int>("matrix fill", 0);
        rtn->set<int>("vector fill", 0);
        rtn->set<int>("assembly loop", 0);
        rtn->set<int>("evaluation", 0);
        rtn->set<int>("evaluation mediator", 0);
        rtn->set("Integration", *ElementIntegral::defaultVerbParams());
        rtn->set("DOF Map", *DOFMapBase::defaultVerbParams());
        first = false;
      }
      return rtn;
    }


private:

  /** */
  void init(const Mesh& mesh, 
    const RefCountPtr<EquationSet>& eqn);

  /** */
  void displayEvaluationResults(
    const EvalContext& context, 
    const EvaluatableExpr* evalExpr, 
    const Array<double>& constantCoeffs, 
    const Array<RefCountPtr<EvalVector> >& vectorCoeffs) const ;

  /** */
  void assemblyLoop(const ComputationType compType,
    RefCountPtr<AssemblyKernelBase> kernel) const ;


  /** */
  void configureMatrix(LinearOperator<double>& A,
    Vector<double>& b) const ;

  /** */
  void configureVector(Vector<double>& b) const ;

  /** */
  void configureMatrixBlock(int br, int bc, 
    LinearOperator<double>& A) const ;

  /** */
  void configureVectorBlock(int br, Vector<double>& b) const ;

  /** */
  Array<Array<int> > findNonzeroBlocks() const ;

  /** */
  IntegrationCellSpecifier whetherToUseCofacets(
    const Array<IntegralGroup>& groups,
    const EvaluatableExpr* ee,
    bool isMaximalCell) const ;

  
  
  /** */
  static int defaultWorkSetSize() {static int rtn=100; return rtn;}

  bool partitionBCs_;
      
  mutable bool matNeedsConfiguration_;
      
  mutable bool matNeedsFinalization_;

  mutable bool vecNeedsConfiguration_;

  Mesh mesh_;

  RefCountPtr<EquationSet> eqn_;

  Array<RefCountPtr<DOFMapBase> > rowMap_;

  Array<RefCountPtr<DOFMapBase> > colMap_;

  Array<RefCountPtr<DiscreteSpace> > externalRowSpace_;

  Array<RefCountPtr<DiscreteSpace> > externalColSpace_;

  Array<RefCountPtr<DiscreteSpace> > privateRowSpace_;

  Array<RefCountPtr<DiscreteSpace> > privateColSpace_;

  Array<RefCountPtr<Set<int> > > bcRows_;

  Array<RefCountPtr<Set<int> > > bcCols_;

  Array<RegionQuadCombo> rqc_;

  Map<ComputationType, Array<EvalContext> > contexts_;

  Array<int> isBCRqc_;

  Map<ComputationType, Array<Array<IntegralGroup> > > groups_;

  Array<RefCountPtr<StdFwkEvalMediator> > mediators_;

  Map<ComputationType, Array<const EvaluatableExpr*> > evalExprs_;

  RefCountPtr<EvalManager> evalMgr_;

  Array<RefCountPtr<Array<int> > > isBCRow_;

  Array<RefCountPtr<Array<int> > > isBCCol_;

  Array<RefCountPtr<std::set<int> > > remoteBCCols_;

  Array<int> lowestRow_;

  Array<int> lowestCol_;

  Array<VectorType<double> > rowVecType_;

  Array<VectorType<double> > colVecType_;

  Map<int, int> testIDToBlockMap_;

  Map<int, int> unkIDToBlockMap_;

  Map<ComputationType, Array<IntegrationCellSpecifier> > rqcRequiresMaximalCofacets_;

  Array<RefCountPtr<PartitionedToMonolithicConverter> > converter_;

};
}
}



#endif
