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

#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceMatrixVectorAssemblyKernel.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;
using std::setw;
using std::endl;
      


void MatrixVectorAssemblyKernel::init(
  const Array<RefCountPtr<DOFMapBase> >& rowMap,
  const Array<RefCountPtr<DOFMapBase> >& colMap,
  LinearOperator<double> A,
  bool partitionBCs)
{
  Tabs tab;
  SUNDANCE_MSG2(verb(), tab << "begin MVAssemblyKernel::init()");

  int numRowBlocks = rowMap.size();
  int numColBlocks = colMap.size();

  for (int br=0; br<numRowBlocks; br++)
    {
      Vector<double> vecBlock; 

      mat_[br].resize(numColBlocks);
      for (int bc=0; bc<numColBlocks; bc++)
        {
          LinearOperator<double> matBlock;
          if (partitionBCs && numRowBlocks==1 && numColBlocks==1)
          {
            matBlock = A;
          }
          else
          {
            matBlock = A.getBlock(br, bc);
          }
          if (matBlock.ptr().get() == 0) continue;
          const Thyra::ZeroLinearOpBase<double>* zp
            = dynamic_cast<const TSFExtended::ZeroLinearOpBase<double>*>(matBlock.ptr().get());
          if (zp) continue;
          mat_[br][bc] 
            = dynamic_cast<TSFExtended::LoadableMatrix<double>* >(matBlock.ptr().get());
          TEST_FOR_EXCEPTION(mat_[br][bc]==0, RuntimeError,
                             "matrix block (" << br << ", " << bc 
                             << ") is not loadable in Assembler::assemble()");
          mat_[br][bc]->zero();
        }
    }
  SUNDANCE_MSG2(verb(), tab << "end MVAssemblyKernel::init()");
}



void MatrixVectorAssemblyKernel::fill(
  bool isBC, 
  const IntegralGroup& group,
  const RefCountPtr<Array<double> >& localValues) 
{
  Tabs tab0;
  SUNDANCE_MSG1(verb(), tab0 << "in MatrixVectorAssemblyKernel::fill()");
  SUNDANCE_MSG1(verb(), tab0 << "filling for integral group " << group.derivs());

  bool useCofacets = group.usesMaximalCofacets();

  if (group.isOneForm())
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb(), tab1 << "group is one form");
    insertLocalVectorBatch(isBC, useCofacets, 
      group.testID(), group.testBlock(), group.mvIndices(),
      *localValues);
  }
  else if (group.isTwoForm())
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb(), tab1 << "group is two form");
    insertLocalMatrixBatch(isBC, useCofacets, 
      group.testID(), group.testBlock(),
      group.unkID(), group.unkBlock(),
      *localValues);
  }
  else
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb(), tab1 << "is zero form -- nothing to do here");
  }

  SUNDANCE_MSG1(verb(), tab0 << "done MatrixVectorAssemblyKernel::fill()");
}
  

void MatrixVectorAssemblyKernel::prepareForWorkSet(
  const Array<Set<int> >& requiredTests,
  const Array<Set<int> >& requiredUnks,
  RefCountPtr<StdFwkEvalMediator> mediator)
{
  Tabs tab0;
  SUNDANCE_MSG1(verb(), tab0 
    << "in MatrixVectorAssemblyKernel::prepareForWorkSet()");

  IntegrationCellSpecifier intCellSpec = mediator->integrationCellSpec();
  SUNDANCE_MSG2(verb(), tab0 
    << "integration cell specifier is " << intCellSpec);
  
  SUNDANCE_MSG2(verb(), tab0 << "building row DOF maps");
  buildLocalDOFMaps(mediator, intCellSpec, requiredTests);

  SUNDANCE_MSG2(verb(), tab0 << "building column DOF maps");
  cmb_.buildLocalDOFMaps(mediator, intCellSpec, requiredUnks);

  SUNDANCE_MSG1(verb(), tab0 << "done MatrixVectorAssemblyKernel::prepareForWorkSet()");
}


void MatrixVectorAssemblyKernel::insertLocalMatrixBatch(
  bool isBCRqc,
  bool useCofacetCells,
  const Array<int>& testID, 
  const Array<int>& testBlock, 
  const Array<int>& unkID,
  const Array<int>& unkBlock,
  const Array<double>& localValues) const
{
  Tabs tab;

  SUNDANCE_MSG1(verb(), tab << "inserting local matrices");

  static Array<int> skipRow;
  static Array<int> rows;
  static Array<int> cols;

  int nCells = rmb().nCells();

  for (unsigned int t=0; t<testID.size(); t++)
  {
    Tabs tab1;
    int br = testBlock[t];

    SUNDANCE_MSG3(verb(), tab1 << "block row = " << br 
      << tab1 << "function ID = "<< testID[t] 
      << " of " << testID.size() << endl
      << tab1 << "is BC eqn = " << isBCRqc << endl
      << tab1 << "num cells = " << nCells << endl
      << tab1 << "using cofacet cells = " << useCofacetCells);

    const RefCountPtr<DOFMapBase>& rowMap = rmb().dofMap(br);
    int lowestLocalRow = rmb().lowestLocalIndex(br);
    int highestRowIndex = lowestLocalRow + rowMap->numLocalDOFs();
    int testChunk = rmb().mapStruct(br, 
      useCofacetCells)->chunkForFuncID(testID[t]);
    int testFuncIndex = rmb().mapStruct(br, 
      useCofacetCells)->indexForFuncID(testID[t]);
    const Array<int>& testDofs = rmb().localDOFs(br, useCofacetCells, testChunk);
    int nTestFuncs = rmb().mapStruct(br, useCofacetCells)->numFuncs(testChunk);
    int nTestNodes = rmb().nNodesInChunk(br, useCofacetCells, testChunk);

    SUNDANCE_MSG3(verb(), tab1 << "lowestLocalRow = " << lowestLocalRow << endl
      << tab1 << "test chunk = " << testChunk << endl
      << tab1 << "func offset into local DOF map = " 
      << testFuncIndex << endl
      << tab1 << "local test dofs = " << testDofs << endl
      << tab1 << "num test funcs in chunk = " << nTestFuncs << endl
      << tab1 << "num test nodes in chunk = " << nTestNodes);
    
    int numRows = nCells * nTestNodes;
    const Array<int>& isBCRow = *(rmb().isBCIndex(br));
    rows.resize(numRows);
    skipRow.resize(numRows);
    int r=0;
    for (int c=0; c<nCells; c++)
    {
      for (int n=0; n<nTestNodes; n++, r++)
      {
        int row = testDofs[(c*nTestFuncs + testFuncIndex)*nTestNodes + n];
        rows[r] = row;
        int localRow = rows[r]-lowestLocalRow;
        skipRow[r] = row < lowestLocalRow || row >= highestRowIndex
          || (isBCRqc && !isBCRow[localRow])
          || (!isBCRqc && isBCRow[localRow]);
      }
    }
    
    for (unsigned int u=0; u<unkID.size(); u++)
    {      
      Tabs tab2;
      int bc = unkBlock[u];
      
      int lowestLocalCol = cmb().lowestLocalIndex(bc);
      int unkChunk = cmb().mapStruct(bc, 
        useCofacetCells)->chunkForFuncID(unkID[u]);
      int unkFuncIndex = cmb().mapStruct(bc, 
        useCofacetCells)->indexForFuncID(unkID[u]);
      const Array<int>& unkDofs = cmb().localDOFs(bc, useCofacetCells, unkChunk);
      int nUnkFuncs = cmb().mapStruct(bc, useCofacetCells)->numFuncs(unkChunk);
      int nUnkNodes = cmb().nNodesInChunk(bc, useCofacetCells, unkChunk);


      SUNDANCE_MSG3(verb(), tab2 << "lowestLocalCol = " 
        << lowestLocalCol << endl
        << tab2 << "block col = " << bc << endl
        << tab2 << "unk ID = "<< unkID[t] 
        << " of " << unkID.size() << endl
        << tab2 << "unk chunk = " << unkChunk << endl
        << tab2 << "func offset into local DOF map = " 
        << unkFuncIndex << endl
        << tab2 << "local unk dofs = " << unkDofs << endl
        << tab2 << "num unk funcs in chunk = " << nUnkFuncs << endl
        << tab2 << "num unk nodes in chunk = " << nUnkNodes);

      cols.resize(nCells*nUnkNodes);

      int j=0;
      for (int c=0; c<nCells; c++)
      {
        for (int n=0; n<nUnkNodes; n++, j++)
        {
          cols[j] = unkDofs[(c*nUnkFuncs + unkFuncIndex)*nUnkNodes + n];
        }
      }
      
      if (verb() >= 3)
      {
        writeLSMs(br, bc, useCofacetCells,
          nTestNodes, nTestFuncs, testFuncIndex, testDofs,
          nUnkNodes, nUnkFuncs, unkFuncIndex, unkDofs, localValues);
      }

      SUNDANCE_MSG2(verb(), tab2 << "calling addToElementBatch()");
      mat_[br][bc]->addToElementBatch(numRows,
        nTestNodes,
        &(rows[0]),
        nUnkNodes,
        &(cols[0]),
        &(localValues[0]),
        &(skipRow[0]));
      SUNDANCE_MSG2(verb(), tab2 << "done calling addToElementBatch()");
    }
  }

  SUNDANCE_MSG1(verb(), tab << "done inserting local matrices");
}                  


void MatrixVectorAssemblyKernel::writeLSMs(
  int blockRow,
  int blockCol,
  bool useCofacetCells,
  int nTestNodes, 
  int nTestFuncs, 
  int testFuncIndex, 
  const Array<int>& testDofs,
  int nUnkNodes, 
  int nUnkFuncs, 
  int unkFuncIndex, 
  const Array<int>& unkDofs,
  const Array<double>& localValues) const
{
  FancyOStream& os = Out::os();

  int nCells = rmb().nCells();
  RefCountPtr<const Array<int> > workSet = rmb().workSet(blockRow, useCofacetCells);

  int lr = 0;

  ios_base::fmtflags oldFlags = os.flags();
  os.setf(ios_base::right);    
  os.setf(ios_base::showpoint);

  for (int c=0; c<nCells; c++)
  {
    Tabs tab3;
    os << tab3 << endl 
       << tab3 << "c="<< c << ", cellLID=" << (*workSet)[c] << endl;

    os << tab3 << "num values per cell = " << localValues.size()/nCells 
       << ", num test nodes=" << nTestNodes << ", num unk nodes="
       << nUnkNodes << endl;
    Array<int> lsmCols(nUnkNodes);
    os << tab3 << setw(17);

    os << "|";
    for (int n=0; n<nUnkNodes; n++)
    {
      int colDof = unkDofs[(c*nUnkFuncs + unkFuncIndex)*nUnkNodes + n];
      lsmCols[n] = colDof;
      os << setw(12) << colDof;
    }
    os << endl << tab3 << "---------------------------------------------------------------" << endl;        

        
    for (int m=0; m<nTestNodes; m++, lr++)
    {
      int globalRow 
        = testDofs[(c*nTestFuncs + testFuncIndex)*nTestNodes+m];
      os << tab3 << setw(6) << m << setw(10) << globalRow << "|";

      for (int n=0; n<nUnkNodes; n++)
      {
        double Amn = localValues[lr*nUnkNodes + n];
        os << setw(12) << setprecision(5) << Amn;
      }
      os << endl;
    }
  }
  os.flags(oldFlags);      
}
