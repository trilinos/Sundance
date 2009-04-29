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

#include "SundanceQuadratureIntegral.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceMultiIndex.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;

extern "C" 
{
  int dgemm_(const char* transA, const char* transB,
             const int* M, const int *N, const int* K,
             const double* alpha, 
             const double* A, const int* ldA,
             const double* B, const int* ldB,
             const double* beta,
             double* C, const int* ldC);
}

static Time& quadratureTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("quadrature"); 
  return *rtn;
}


QuadratureIntegral::QuadratureIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim, 
  const CellType& cellType,
  const QuadratureFamily& quad,
  const ParameterList& verbParams)
  : QuadratureIntegralBase(spatialDim, maxCellType, dim, cellType, quad , verbParams),
    W_(),
    useSumFirstMethod_(true)
{
  Tabs tab0;
  int verb = verbLevel("setup");
  SUNDANCE_MSG1(verb, tab0 << "QuadratureIntegral ctor for zero-form");
  Tabs tab1;
  SUNDANCE_MSG1(verb, tab1 << "cellDim=" << dim 
    << ", spatialDim=" << spatialDim);
  SUNDANCE_MSG1(verb, tab1 << "quadrature family=" << quad);  


  W_.resize(nFacetCases());

  for (int fc=0; fc<nFacetCases(); fc++)
    {
      Tabs tab2;
      SUNDANCE_MSG1(verb, tab2 << "facet case=" << fc
        << " of " << nFacetCases());        
      Tabs tab3;
      /* create the quad points and weights */
      Array<double> quadWeights;
      Array<Point> quadPts;
      if (nFacetCases()==1) quad.getPoints(cellType, quadPts, quadWeights);
      else quad.getFacetPoints(maxCellType, dim, fc, quadPts, quadWeights);
      nQuad_ = quadPts.size();
      
      W_[fc].resize(nQuad());
      
      SUNDANCE_MSG1(verb, tab3 << "num quad pts" << nQuad());
      
      SUNDANCE_MSG1(verb, tab3 << "quad weights" << quadWeights);
      
      for (int q=0; q<nQuad(); q++)
        {
          W_[fc][q] = quadWeights[q];
        }
    }    
}


QuadratureIntegral::QuadratureIntegral(int spatialDim,
                                       const CellType& maxCellType,
                                       int dim, 
                                       const CellType& cellType,
                                       const BasisFamily& testBasis,
                                       int alpha,
                                       int testDerivOrder,
                                       const QuadratureFamily& quad,
				       const ParameterList& verbParams)
  : QuadratureIntegralBase(spatialDim, maxCellType, dim, cellType, 
			   testBasis, alpha, testDerivOrder, quad , verbParams),
    W_(),
    useSumFirstMethod_(true)
{
  Tabs tab0;
  int verb = verbLevel("setup");
  SUNDANCE_MSG1(verb, tab0 << "QuadratureIntegral ctor for one-form");
  Tabs tab1;
  SUNDANCE_MSG1(verb, tab1 << "cellDim=" << dim 
    << ", spatialDim=" << spatialDim);
  SUNDANCE_MSG1(verb, tab1 << "quadrature family=" << quad);  
  SUNDANCE_MSG1(verb,
               tab1 
               << "******** computing basis functions on quad pts *******"
               << std::endl << tab0 << "test basis=" 
               << testBasis 
               << std::endl << tab0 << "cell type=" << cellType
               << std::endl << tab0 << "differentiation order="
               << testDerivOrder);
  SUNDANCE_MSG1(verb,
               tab1 << "num test derivative components=" 
               << nRefDerivTest());

  TEST_FOR_EXCEPTION(testDerivOrder < 0 || testDerivOrder > 1,
                     InternalError,
                     "Test function derivative order=" << testDerivOrder
                     << " must be 0 or 1");
  

  W_.resize(nFacetCases());

  CellType evalCellType = cellType;
  if (nFacetCases() > 1) evalCellType = maxCellType;

  SUNDANCE_MSG1(verb, "evaluation cell type is " << evalCellType);

  for (int fc=0; fc<nFacetCases(); fc++)
    {
      Tabs tab2;

      SUNDANCE_MSG1(verb, tab2 << "facet case=" << fc
        << " of " << nFacetCases());        

      /* create the quad points and weights */
      Array<double> quadWeights;
      Array<Point> quadPts;
      if (nFacetCases()==1) quad.getPoints(cellType, quadPts, quadWeights);
      else quad.getFacetPoints(maxCellType, dim, fc, quadPts, quadWeights);
      nQuad_ = quadPts.size();
      
      W_[fc].resize(nQuad() * nRefDerivTest() * nNodesTest());

      SUNDANCE_MSG1(verb, tab2 << "num quad pts " << nQuad());
      SUNDANCE_MSG1(verb, tab2 << "num nodes for test function " << nNodesTest());

      Array<Array<Array<Array<double> > > > testBasisVals(nRefDerivTest());


      for (int r=0; r<nRefDerivTest(); r++)
        {
          Tabs tab3;
          SUNDANCE_MSG1(verb, tab3 
            << "evaluating basis functions for ref deriv direction " << r);
          MultiIndex mi;
          testBasisVals[r].resize(testBasis.dim());
          if (testDerivOrder==1) mi[r] = 1;
          testBasis.ptr()->refEval(maxCellType, evalCellType, quadPts, mi, 
                                   testBasisVals[r]);
          for (unsigned int c=0; c<testBasisVals[r].size(); c++)
          {
            Tabs tab4;
            SUNDANCE_MSG1(verb, tab4<< "vector component c=" << c 
              << " of " << testBasisVals[r].size());
            for (unsigned int q=0; q<testBasisVals[r][c].size(); q++)
            {
              Tabs tab5;
              SUNDANCE_MSG1(verb, tab5<< "quad point=" << q 
              << " basis values= " << testBasisVals[r][c][q]);
            }
          }
        }

      
      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "all basis values " << testBasisVals);

      int vecComp = 0;
      for (int q=0; q<nQuad(); q++)
        {
          for (int t=0; t<nRefDerivTest(); t++)
            {
              for (int nt=0; nt<nNodesTest(); nt++)
                {
                  wValue(fc, q, t, nt) 
                    = chop(quadWeights[q] * testBasisVals[t][vecComp][q][nt]) ;
                }
            }
        }    
  
      addFlops(2*nQuad()*nRefDerivTest()*nNodesTest());
    }
  if (verbosity() > VerbMedium)
    {
      print(std::cerr);
    }
}




QuadratureIntegral::QuadratureIntegral(int spatialDim,
                                       const CellType& maxCellType,
                                       int dim,
                                       const CellType& cellType,
                                       const BasisFamily& testBasis,
                                       int alpha,
                                       int testDerivOrder,
                                       const BasisFamily& unkBasis,
                                       int beta,
                                       int unkDerivOrder,
                                       const QuadratureFamily& quad,
				       const ParameterList& verbParams)
  : QuadratureIntegralBase(spatialDim, maxCellType, dim, cellType, 
			   testBasis, alpha, testDerivOrder, 
			   unkBasis, beta, unkDerivOrder, quad , verbParams), 
    W_(),
    useSumFirstMethod_(true)
{
  Tabs tab0;
  verbosity() = classVerbosity();
  SUNDANCE_OUT(this->verbosity() > VerbSilent, 
               tab0 << " ************* computing basis func products on quad pts ***************" 
               << std::endl << tab0 << "test basis=" 
               << testBasis 
               << std::endl << tab0 << "unk basis=" << unkBasis
               << std::endl << tab0 << "cell type=" << cellType
               << std::endl << tab0 <<"differentiation order t=" 
               << testDerivOrder << ", u=" << unkDerivOrder);
  SUNDANCE_OUT(this->verbosity() > VerbMedium, 
               tab0 << "num test derivative components=" 
               << nRefDerivTest());

  TEST_FOR_EXCEPTION(testDerivOrder < 0 || testDerivOrder > 1,
                     InternalError,
                     "Test function derivative order=" << testDerivOrder
                     << " must be 0 or 1");
  
  TEST_FOR_EXCEPTION(unkDerivOrder < 0 || unkDerivOrder > 1,
                     InternalError,
                     "Unknown function derivative order=" << unkDerivOrder
                     << " must be 0 or 1");

  

  W_.resize(nFacetCases());

  CellType evalCellType = cellType;
  if (nFacetCases() > 1) evalCellType = maxCellType;

  for (int fc=0; fc<nFacetCases(); fc++)
    {
      /* get the quad pts and weights */
      Array<double> quadWeights;
      Array<Point> quadPts;
      if (nFacetCases()==1) quad.getPoints(cellType, quadPts, quadWeights);
      else quad.getFacetPoints(maxCellType, dim, fc, quadPts, quadWeights);
      nQuad_ = quadPts.size();

      W_[fc].resize(nQuad() * nRefDerivTest() * nNodesTest()  
                * nRefDerivUnk() * nNodesUnk());

      SUNDANCE_OUT(this->verbosity() > VerbLow, 
                   tab0 << "num quad pts" << nQuad());

      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "quad pts" << quadPts);

      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "quad weights" << quadWeights);



      /* compute the basis functions */
      Array<Array<Array<Array<double> > > > testBasisVals(nRefDerivTest());
      Array<Array<Array<Array<double> > > > unkBasisVals(nRefDerivUnk());


      for (int r=0; r<nRefDerivTest(); r++)
        {
          testBasisVals[r].resize(testBasis.dim());
          MultiIndex mi;
          if (testDerivOrder==1) mi[r] = 1;
          testBasis.ptr()->refEval(maxCellType, evalCellType, quadPts, mi, 
                                   testBasisVals[r]);
        }

      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "test basis values" << testBasisVals);

      for (int r=0; r<nRefDerivUnk(); r++)
        {
          unkBasisVals[r].resize(unkBasis.dim());
          MultiIndex mi;
          if (unkDerivOrder==1) mi[r] = 1;
          unkBasis.ptr()->refEval(maxCellType, evalCellType, 
                                  quadPts, mi, unkBasisVals[r]);
        }


      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "unk basis values" << unkBasisVals);


      int vecComp = 0;
      /* form the products of basis functions at each quad pt */
      for (int q=0; q<nQuad(); q++)
        {
          for (int t=0; t<nRefDerivTest(); t++)
            {
              for (int nt=0; nt<nNodesTest(); nt++)
                {
                  for (int u=0; u<nRefDerivUnk(); u++)
                    {
                      for (int nu=0; nu<nNodesUnk(); nu++)
                        {
                          wValue(fc, q, t, nt, u, nu)
                            = chop(quadWeights[q] * testBasisVals[t][vecComp][q][nt] 
                                   * unkBasisVals[u][vecComp][q][nu]);
                        }
                    }
                }
            }
        }

      addFlops(3*nQuad()*nRefDerivTest()*nNodesTest()*nRefDerivUnk()*nNodesUnk()
               + W_[fc].size());
      for (unsigned int i=0; i<W_[fc].size(); i++) W_[fc][i] = chop(W_[fc][i]);
    
    }
  if (verbosity() > VerbMedium)
    {
      print(std::cerr);
    }

}


void QuadratureIntegral::transformZeroForm(const CellJacobianBatch& JTrans,  
                                           const CellJacobianBatch& JVol,
                                           const Array<int>& isLocalFlag,
                                           const Array<int>& facetIndex,
                                           const double* const coeff,
                                           RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(quadratureTimer());
  Tabs tabs;
  TEST_FOR_EXCEPTION(order() != 0, InternalError,
                     "QuadratureIntegral::transformZeroForm() called "
                     "for form of order " << order());

  TEST_FOR_EXCEPTION( (int) isLocalFlag.size() != 0 
                      && (int) isLocalFlag.size() != JVol.numCells(),
                      RuntimeError,
                      "mismatch between isLocalFlag.size()=" << isLocalFlag.size()
                      << " and JVol.numCells()=" << JVol.numCells());

  bool checkLocalFlag = (int) isLocalFlag.size() != 0;

  
  SUNDANCE_VERB_MEDIUM(tabs << "doing zero form by quadrature");
  double& a = (*A)[0];
  SUNDANCE_VERB_EXTREME(tabs << "input A = " << *A);
  double* coeffPtr = (double*) coeff;

  int flops = 0;
  if (nFacetCases()==1)
    {
      const Array<double>& w = W_[0];
      for (int c=0; c<JVol.numCells(); c++)
        {
          if (checkLocalFlag && !isLocalFlag[c]) 
            {
              coeffPtr += nQuad();
              continue;
            }
          double detJ = fabs(JVol.detJ()[c]);
          for (int q=0; q<nQuad(); q++, coeffPtr++)
            {
              a += w[q]*(*coeffPtr)*detJ;
            }
          flops += 3*nQuad();
        }
    }
  else
    {
      for (int c=0; c<JVol.numCells(); c++)
        {
          if (checkLocalFlag && !isLocalFlag[c]) 
            {
              coeffPtr += nQuad();
              continue;
            }
          double detJ = fabs(JVol.detJ()[c]);
          const Array<double>& w = W_[facetIndex[c]];
          for (int q=0; q<nQuad(); q++, coeffPtr++)
            {
              a += w[q]*(*coeffPtr)*detJ;
            }
          flops += 3*nQuad();
        }
    }
  SUNDANCE_VERB_EXTREME(tabs << "output A = " << *A);

  addFlops(flops);
}



void QuadratureIntegral::transformOneForm(const CellJacobianBatch& JTrans,  
                                          const CellJacobianBatch& JVol,
                                          const Array<int>& facetIndex,
                                          const double* const coeff,
                                          RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(quadratureTimer());
  Tabs tabs;
  TEST_FOR_EXCEPTION(order() != 1, InternalError,
                     "QuadratureIntegral::transformOneForm() called for form "
                     "of order " << order());
  SUNDANCE_VERB_MEDIUM(tabs << "doing one form by quadrature");
  int flops = 0;

  /* If the derivative order is zero, the only thing to be done 
   * is to multiply by the cell's Jacobian determinant and sum over the
   * quad points */
  if (testDerivOrder() == 0)
    {
      double* aPtr = &((*A)[0]);
      SUNDANCE_VERB_EXTREME(tabs << "input A = " << *A);
      double* coeffPtr = (double*) coeff;
      int offset = 0 ;
      const Array<double>& w = W_[0];

      for (int c=0; c<JVol.numCells(); c++, offset+=nNodes())
        {
          Tabs tab2;
          double detJ = fabs(JVol.detJ()[c]);
          SUNDANCE_VERB_EXTREME(tab2 << "c=" << c << " detJ=" << detJ);
          for (int q=0; q<nQuad(); q++, coeffPtr++)
            {
              Tabs tab3;
              double f = (*coeffPtr)*detJ;
              SUNDANCE_VERB_EXTREME(tab3 << "q=" << q << " coeff=" << 
                                    *coeffPtr << " f=" << f);
              for (int n=0; n<nNodes(); n++) 
                {
                  Tabs tab4;
                  SUNDANCE_VERB_EXTREME(tab4 << "n=" << n << " w=" <<  
                                        w[n + nNodes()*q]);
                  aPtr[offset+n] += f*w[n + nNodes()*q];
                }
            }
        }
      
      SUNDANCE_VERB_EXTREME(tabs << "output A = " << *A);
      addFlops( JVol.numCells() * (1 + nQuad() * (1 + 2*nNodes())) );
    }
  else
    {
      /* If the derivative order is nonzero, then we have to do a transformation. 
       * If we're also on a cell of dimension lower than maximal, we need to refer
       * to the facet index of the facet being integrated. */

      createOneFormTransformationMatrix(JTrans, JVol);

      SUNDANCE_OUT(this->verbosity() > VerbMedium, 
                   Tabs() << "transformation matrix=" << G(alpha()));

      double* GPtr = &(G(alpha())[0]);      

      if (useSumFirstMethod())
        {
          transformSummingFirst(JVol.numCells(), facetIndex, GPtr, coeff, A);
        }
      else
        {
          transformSummingLast(JVol.numCells(), facetIndex, GPtr, coeff, A);
        }
    }
  addFlops(flops);
}


void QuadratureIntegral::transformTwoForm(const CellJacobianBatch& JTrans,
                                          const CellJacobianBatch& JVol,
                                          const Array<int>& facetIndex,
                                          const double* const coeff,
                                          RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(quadratureTimer());
  Tabs tabs;
  TEST_FOR_EXCEPTION(order() != 2, InternalError,
                     "QuadratureIntegral::transformTwoForm() called for form "
                     "of order " << order());
  SUNDANCE_VERB_MEDIUM(tabs << "doing one form by quadrature");

  /* If the derivative orders are zero, the only thing to be done 
   * is to multiply by the cell's Jacobian determinant and sum over the
   * quad points */
  if (testDerivOrder() == 0 && unkDerivOrder() == 0)
    {
      double* aPtr = &((*A)[0]);
      double* coeffPtr = (double*) coeff;
      int offset = 0 ;

      const Array<double>& w = W_[0];
      for (int c=0; c<JVol.numCells(); c++, offset+=nNodes())
        {
          double detJ = fabs(JVol.detJ()[c]);
          for (int q=0; q<nQuad(); q++, coeffPtr++)
            {
              double f = (*coeffPtr)*detJ;
              for (int n=0; n<nNodes(); n++) 
                {
                  aPtr[offset+n] += f*w[n + nNodes()*q];
                }
            }
        }

      addFlops( JVol.numCells() * (1 + nQuad() * (1 + 2*nNodes())) );
    }
  else
    {
      createTwoFormTransformationMatrix(JTrans, JVol);
      double* GPtr;

      if (testDerivOrder() == 0)
        {
          GPtr = &(G(beta())[0]);
          SUNDANCE_OUT(this->verbosity() > VerbMedium, 
                       Tabs() << "transformation matrix=" << G(beta()));
        }
      else if (unkDerivOrder() == 0)
        {
          GPtr = &(G(alpha())[0]);
          SUNDANCE_OUT(this->verbosity() > VerbMedium, 
                       Tabs() << "transformation matrix=" << G(alpha()));
        }
      else
        {
          GPtr = &(G(alpha(), beta())[0]);
          SUNDANCE_OUT(this->verbosity() > VerbMedium, 
                       Tabs() << "transformation matrix=" 
                       << G(alpha(),beta()));
        }
        
      
      if (useSumFirstMethod())
        {
          transformSummingFirst(JTrans.numCells(), facetIndex, GPtr, coeff, A);
        }
      else
        {
          transformSummingLast(JTrans.numCells(), facetIndex, GPtr, coeff, A);
        }
    }
}

void QuadratureIntegral
::transformSummingFirst(int nCells,
                        const Array<int>& facetIndex,
                        const double* const GPtr,
                        const double* const coeff,
                        RefCountPtr<Array<double> >& A) const
{
  double* aPtr = &((*A)[0]);
  double* coeffPtr = (double*) coeff;
  
  int transSize = 0; 
  if (order()==2)
    {
      transSize = nRefDerivTest() * nRefDerivUnk();
    }
  else
    {
      transSize = nRefDerivTest();
    }

  if (false)
    {
      std::cerr << "nCells = " << nCells << std::endl;
      std::cerr << "nNodes = " << nNodes() << std::endl;
      std::cerr << "nQuad = " << nQuad() << std::endl;
      std::cerr << "transSize = " << transSize << std::endl;
    }
  /* The sum workspace is used to store the sum of untransformed quantities */
  static Array<double> sumWorkspace;

  int swSize = transSize * nNodes();
  sumWorkspace.resize(swSize);
  
  
  /*
   * The number of operations for the sum-first method is 
   * 
   * Adds: (N_c * nNodes * transSize) * (N_q + 1) 
   * Multiplies: same as number of adds
   * Total: 2*(N_c * nNodes * transSize) * (N_q + 1) 
   */
  
  for (int c=0; c<nCells; c++)
    {
      /* sum untransformed basis combinations over quad points */
      for (int i=0; i<swSize; i++) sumWorkspace[i]=0.0;
      int fc = 0;
      if (nFacetCases() > 1) fc = facetIndex[c];

      const Array<double>& w = W_[fc];
      for (int q=0; q<nQuad(); q++, coeffPtr++)
        {
          double f = (*coeffPtr);
          for (int n=0; n<swSize; n++) 
            {
              sumWorkspace[n] += f*w[n + q*swSize];
            }
        }
      /* transform the sum */
      const double * const gCell = &(GPtr[transSize*c]);
      double* aCell = aPtr + nNodes()*c;
      for (int i=0; i<nNodes(); i++)
        {
          for (int j=0; j<transSize; j++)
            {
              aCell[i] += sumWorkspace[nNodes()*j + i] * gCell[j];
            }
        }
    }
  
  int flops = 2*(nCells * nNodes() * transSize) * (nQuad() + 1) ;
  addFlops(flops);
}

void QuadratureIntegral
::transformSummingLast(int nCells,
                       const Array<int>& facetIndex,
                       const double* const GPtr,
                       const double* const coeff,
                       RefCountPtr<Array<double> >& A) const
{
  double* aPtr = &((*A)[0]);
  int transSize = 0; 
  
  if (order()==2)
    {
      transSize = nRefDerivTest() * nRefDerivUnk();
    }
  else
    {
      transSize = nRefDerivTest();
    }

  /* This workspace is used to store the jacobian values scaled by the coeff
   * at that quad point */
  static Array<double> jWorkspace;
  jWorkspace.resize(transSize);


  /*
   * The number of operations for the sum-last method is 
   * 
   * Adds: N_c * N_q * transSize * nNodes
   * Multiplies: numAdds + N_c * N_q * transSize
   * Total: N_c * N_q * transSize * (1 +  2*nNodes)
   */

  for (int c=0; c<nCells; c++)
    {
      const double* const gCell = &(GPtr[transSize*c]);
      double* aCell = aPtr + nNodes()*c;
      int fc = 0;
      if (nFacetCases() > 1) fc = facetIndex[c];
      const Array<double>& w = W_[fc];

      for (int q=0; q<nQuad(); q++)
        {
          double f = coeff[c*nQuad() + q];
          for (int t=0; t<transSize; t++) jWorkspace[t]=f*gCell[t];

          for (int n=0; n<nNodes(); n++)
            {
              for (int t=0; t<transSize; t++)
                {
                  aCell[n] += jWorkspace[t]*w[n + nNodes()*(t + transSize*q)];
                }
            }
        }
    }
  int flops = nCells * nQuad() * transSize * (1 + 2*nNodes()) ;
  addFlops(flops);
}

