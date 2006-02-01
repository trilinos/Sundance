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
#include "TSFIncrementallyConfigurableMatrixFactory.hpp"
#include "TSFCollectivelyConfigurableMatrixFactory.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  namespace Internal
  {
    using namespace Teuchos;

    typedef std::set<int> ColSetType;

    /** 
     * 
     */
    class Assembler : public TSFExtended::ObjectWithVerbosity<Assembler>
    {
    public:
      /** */
      Assembler(const Mesh& mesh, 
                const RefCountPtr<EquationSet>& eqn,
                const Array<VectorType<double> >& rowVectorType,
                const Array<VectorType<double> >& colVectorType,
                const VerbositySetting& verb = classVerbosity());
      /** */
      Assembler(const Mesh& mesh, 
                const RefCountPtr<EquationSet>& eqn,
                const VerbositySetting& verb = classVerbosity());
      
      /** */
      const Array<RefCountPtr<DOFMapBase> >& rowMap() const 
      {return rowMap_;}

      /** */
      const Array<RefCountPtr<DOFMapBase> >& colMap() const 
      {return colMap_;}

      /** */
      const Array<RefCountPtr<DiscreteSpace> >& solutionSpace() const 
      {return colSpace_;}

      /** */
      const Array<RefCountPtr<DiscreteSpace> >& rowSpace() const 
      {return rowSpace_;}

      /** */
      VectorSpace<double> solnVecSpace() const ;

      /** */
      VectorSpace<double> rowVecSpace() const ;

      /** */
      const Array<RefCountPtr<Set<int> > >& bcRows() {return bcRows_;}

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
      static unsigned int& workSetSize() 
      {static unsigned int rtn = defaultWorkSetSize(); return rtn;}
      
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
      static int& numAssembleCalls() {static int rtn=0; return rtn;}

      /** */
      static bool& matrixEliminatesRepeatedCols() {static bool x = false; return x;}
      
    private:

      /** */
      void init(const Mesh& mesh, 
                const RefCountPtr<EquationSet>& eqn);

      /** */
      void insertLocalMatrixBatch(int cellDim, 
                                  const Array<int>& workSet, 
                                  bool isBCRqc,
                                  const Array<Array<Array<int> > >& testIndices,
                                  const Array<Array<Array<int> > >& unkIndices,
                                  const Array<Array<int> >& nTestNodes, 
                                  const Array<Array<int> >& nUnkNodes,
                                  const Array<int>& testID, 
                                  const Array<int>& testBlock, 
                                  const Array<int>& unkID,
                                  const Array<int>& unkBlock,
                                  const Array<double>& localValues, 
                                  Array<Array<LoadableMatrix<double>* > >& mat) const ;

      /** */
      void insertLocalVectorBatch(int cellDim, 
                                  const Array<int>& workSet, 
                                  bool isBCRqc,
                                  const Array<Array<Array<int> > >& testIndices,
                                  const Array<Array<int> >& nTestNodes, 
                                  const Array<int>& testID, 
                                  const Array<int>& testBlock, 
                                  const Array<double>& localValues, 
                                  Array<TSFExtended::LoadableVector<double>* >& vec) const ;

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

      Array<Array<int> > findNonzeroBlocks() const ;

      /** */
      static int defaultWorkSetSize() {return 100;}
      
      mutable bool matNeedsConfiguration_;
      
      mutable bool matNeedsFinalization_;

      mutable bool vecNeedsConfiguration_;

      Mesh mesh_;

      RefCountPtr<EquationSet> eqn_;

      Array<RefCountPtr<DOFMapBase> > rowMap_;

      Array<RefCountPtr<DOFMapBase> > colMap_;

      Array<RefCountPtr<DiscreteSpace> > rowSpace_;

      Array<RefCountPtr<DiscreteSpace> > colSpace_;

      Array<RefCountPtr<Set<int> > > bcRows_;

      Array<RegionQuadCombo> rqc_;

      Map<ComputationType, Array<EvalContext> > contexts_;

      Array<int> isBCRqc_;

      Map<ComputationType, Array<Array<IntegralGroup> > > groups_;

      Array<RefCountPtr<StdFwkEvalMediator> > mediators_;

      Map<ComputationType, Array<const EvaluatableExpr*> > evalExprs_;

      RefCountPtr<EvalManager> evalMgr_;

      Array<RefCountPtr<Array<int> > > isBCRow_;

      Array<int> lowestRow_;

      Array<VectorType<double> > rowVecType_;

      Array<VectorType<double> > colVecType_;

      Map<int, int> testIDToBlockMap_;

      Map<int, int> unkIDToBlockMap_;

    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
