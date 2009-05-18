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

#include "SundanceRefIntegral.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceGaussianQuadratureType.hpp"
#include "SundanceQuadratureType.hpp"
#include "SundanceMultiIndex.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;
using std::ios_base;
using std::setw;
using std::endl;

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

static Time& refIntegrationTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("ref integration"); 
  return *rtn;
}



RefIntegral::RefIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim, 
  const CellType& cellType,
  int verb)
  : ElementIntegral(spatialDim, maxCellType, dim, cellType, verb), W_()
{
  Tabs tab0(0);

  SUNDANCE_MSG1(setupVerb(),
    tab0 << "************* creating reference 0-form integrals ********");
  if (setupVerb()) describe(Out::os());
  
  /* we need to sum the quadrature weights 
     to compute the volume of the reference cell */
  QuadratureFamily quad = new GaussianQuadrature(2);
  Array<Point> quadPts;
  Array<double> quadWeights;

  W_.resize(1);
  W_[0].resize(1);

  quad.getPoints(cellType, quadPts, quadWeights);  

  for (unsigned int q=0; q<quadWeights.size(); q++) W_[0][0] += quadWeights[q];
}

RefIntegral::RefIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim, 
  const CellType& cellType,
  const BasisFamily& testBasis,
  int alpha,
  int testDerivOrder,
  int verb)
  : ElementIntegral(spatialDim, maxCellType, dim, cellType, 
    testBasis, alpha, testDerivOrder, verb), W_()
{
  Tabs tab0(0);
  SUNDANCE_MSG1(setupVerb(),
    tab0 << "************* creating reference 1-form integrals ********");
  if (setupVerb()) describe(Out::os());
  assertLinearForm();

  W_.resize(nFacetCases());

  /* Determine the quadrature order needed for exact integrations */
  QuadratureType qType = new GaussianQuadratureType();
  int reqOrder = qType.findValidOrder(cellType, 
    std::max(1, testBasis.order()));

  SUNDANCE_MSG2(setupVerb(),
    tab0 << "using quadrature order=" << reqOrder);
   
  /* Create a quadrature family of the required order */
  QuadratureFamily quad = qType.createQuadFamily(reqOrder);
  
  /* We now loop over the different evaluation cases, integrating the
   * basis functions for each. Because this is a reference integral,
   * we can actually do the untransformed integrals here. */
  for (int fc=0; fc<nFacetCases(); fc++)
  {
    Tabs tab1;
    SUNDANCE_MSG2(setupVerb(),
      tab1 << "evaluation case=" << fc << " of " << nFacetCases());
    /* initialize size of untransformed integral results array */
    W_[fc].resize(nRefDerivTest() * nNodesTest());

    /* initialize values of integrals to zero */
    for (unsigned int i=0; i<W_[fc].size(); i++) W_[fc][i]=0.0;

    Array<Array<Array<Array<double> > > > testBasisVals(nRefDerivTest());
  
    /* get quadrature points */
    Array<Point> quadPts;
    Array<double> quadWeights;
    getQuad(quad, fc, quadPts, quadWeights);

    int nQuad = quadPts.size();

    /* compute the basis functions */
    for (int r=0; r<nRefDerivTest(); r++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(setupVerb(), tab2 << "evaluating basis derivative " 
        << r << " of " << nRefDerivTest());

      testBasisVals[r].resize(testBasis.dim());
      MultiIndex mi;
      if (testDerivOrder==1) mi[r] = 1;
      testBasis.refEval(evalCellType(), quadPts, mi, 
        testBasisVals[r], setupVerb());
    }

    /* do the quadrature */
    SUNDANCE_MSG2(setupVerb(), tab1 << "doing quadrature");
    int vecComp = 0;
    for (int q=0; q<nQuad; q++)
    {
      for (int t=0; t<nRefDerivTest(); t++)
      {
        for (int nt=0; nt<nNodesTest(); nt++)
        {
          value(fc, t, nt) 
            += chop(quadWeights[q] * testBasisVals[t][vecComp][q][nt]) ;
        }
      }
    }    

    for (unsigned int i=0; i<W_[fc].size(); i++) W_[fc][i] = chop(W_[fc][i]);

    addFlops(3*nQuad*nRefDerivTest()*nNodesTest() + W_[fc].size());
  }

  /* print the result */
  SUNDANCE_MSG4(setupVerb(), tab0 << "--------------------------------------");
  SUNDANCE_MSG4(setupVerb(), tab0 << "reference linear form integral results");
  if (setupVerb() >= 4)
  {
    for (int fc=0; fc<nFacetCases(); fc++)
    {
      Tabs tab1;
      SUNDANCE_MSG4(setupVerb(), tab1 << "------ evaluation case " << fc << " of "
        << nFacetCases() << "-------");
      
      for (int r=0; r<nRefDerivTest(); r++)
      {
        Tabs tab2;

        MultiIndex mi;
        if (testDerivOrder==1) mi[r] = 1;
        SUNDANCE_MSG1(setupVerb(), tab2 << "multiindex=" << mi);

        ios_base::fmtflags oldFlags = Out::os().flags();
        Out::os().setf(ios_base::right);    
        Out::os().setf(ios_base::showpoint);
        for (int nt=0; nt<nNodesTest(); nt++)
        {
          Tabs tab3;
          Out::os() << tab3 << setw(10) << nt 
                    << setw(12) << std::setprecision(5) << value(fc, r, nt) 
                    << endl;
        }
        Out::os().flags(oldFlags);
      }
    }
  }

  SUNDANCE_MSG1(setupVerb(), tab0 << "done reference linear form ctor");
}




RefIntegral::RefIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim,
  const CellType& cellType,
  const BasisFamily& testBasis,
  int alpha,
  int testDerivOrder,
  const BasisFamily& unkBasis,
  int beta,
  int unkDerivOrder, 
  int verb)
  : ElementIntegral(spatialDim, maxCellType,  dim, cellType,
    testBasis, alpha, testDerivOrder, 
    unkBasis, beta, unkDerivOrder, verb), W_()

{
  Tabs tab0(0);
  SUNDANCE_MSG1(setupVerb(),
    tab0 << "************* creating reference 2-form integrals ********");
  if (setupVerb()) describe(Out::os());

  assertBilinearForm();

  W_.resize(nFacetCases());

  QuadratureType qType = new GaussianQuadratureType();
  int reqOrder = qType.findValidOrder(cellType, 
    std::max(1, unkBasis.order() + testBasis.order()));

  SUNDANCE_MSG2(setupVerb(),
    tab0 << "using quadrature order=" << reqOrder);
  QuadratureFamily quad = qType.createQuadFamily(reqOrder);


  SUNDANCE_MSG2(setupVerb(),
    tab0 << "processing evaluation cases");

  for (int fc=0; fc<nFacetCases(); fc++)
  {
    Tabs tab1;
    SUNDANCE_MSG1(setupVerb(), tab1 << "------ evaluation case " << fc << " of "
      << nFacetCases() << "-------");
    
    W_[fc].resize(nRefDerivTest() * nNodesTest()  * nRefDerivUnk() * nNodesUnk());
    for (unsigned int i=0; i<W_[fc].size(); i++) W_[fc][i]=0.0;

    Array<Array<Array<Array<double> > > > testBasisVals(nRefDerivTest());
    Array<Array<Array<Array<double> > > > unkBasisVals(nRefDerivUnk());
        
    Array<Point> quadPts;
    Array<double> quadWeights;
    getQuad(quad, fc, quadPts, quadWeights);
    int nQuad = quadPts.size();

    for (int r=0; r<nRefDerivTest(); r++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(setupVerb(), tab2 
        << "evaluating test function basis derivative " 
        << r << " of " << nRefDerivTest());
      testBasisVals[r].resize(testBasis.dim());
      MultiIndex mi;
      if (testDerivOrder==1) mi[r] = 1;
      testBasis.refEval(evalCellType(), quadPts, mi, 
        testBasisVals[r], setupVerb());
    }

    for (int r=0; r<nRefDerivUnk(); r++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(setupVerb(), tab2 
        << "evaluating unknown function basis derivative " 
        << r << " of " << nRefDerivUnk());
      unkBasisVals[r].resize(unkBasis.dim());
      MultiIndex mi;
      if (unkDerivOrder==1) mi[r] = 1;
      unkBasis.refEval(evalCellType(), 
        quadPts, mi, unkBasisVals[r], setupVerb());
    }

    SUNDANCE_MSG2(setupVerb(), tab1 << "doing quadrature...");
    int vecComp = 0;
    for (int q=0; q<nQuad; q++)
    {
      for (int t=0; t<nRefDerivTest(); t++)
      {
        for (int nt=0; nt<nNodesTest(); nt++)
        {
          for (int u=0; u<nRefDerivUnk(); u++)
          {
            for (int nu=0; nu<nNodesUnk(); nu++)
            {
              value(fc, t, nt, u, nu) 
                += chop(quadWeights[q] * testBasisVals[t][vecComp][q][nt] 
                  * unkBasisVals[u][vecComp][q][nu]);
            }
          }
        }
      }
    }
    SUNDANCE_MSG2(setupVerb(), tab1 << "...done");
    addFlops(4*nQuad*nRefDerivTest()*nNodesTest()*nRefDerivUnk()*nNodesUnk()
      + W_[fc].size());
    for (unsigned int i=0; i<W_[fc].size(); i++) W_[fc][i] = chop(W_[fc][i]);
  }

  SUNDANCE_MSG1(setupVerb(), tab0 
    << "----------------------------------------");
  SUNDANCE_MSG4(setupVerb(), tab0 
    << "reference bilinear form integral results");
  if (setupVerb() >= 4)
  {
    for (int fc=0; fc<nFacetCases(); fc++)
    {
      Tabs tab1;
      SUNDANCE_MSG4(setupVerb(), tab1 << "evaluation case " << fc << " of "
        << nFacetCases());
      
      for (int rt=0; rt<nRefDerivTest(); rt++)
      {
        for (int ru=0; ru<nRefDerivUnk(); ru++)
        {
          Tabs tab2;
          MultiIndex miTest;
          if (testDerivOrder==1) miTest[rt] = 1;
          MultiIndex miUnk;
          if (unkDerivOrder==1) miUnk[ru] = 1;
          SUNDANCE_MSG1(setupVerb(), tab2 << "test multiindex=" << miTest
            << " unk multiindex=" << miUnk);
          
          ios_base::fmtflags oldFlags = Out::os().flags();
          Out::os().setf(ios_base::right);    
          Out::os().setf(ios_base::showpoint);
          for (int nt=0; nt<nNodesTest(); nt++)
          {
            Tabs tab3;
            Out::os() << tab3 << setw(10) << nt;
            for (int nu=0; nu<nNodesUnk(); nu++)
            {
              Out::os() << setw(12) << std::setprecision(5)
                        << value(fc, rt, nt, ru, nu) ;
            }
            Out::os() << endl;
          }
          Out::os().flags(oldFlags);
        }
      }
    }
  }

  SUNDANCE_MSG1(setupVerb(), tab0 << "done reference bilinear form ctor");
}




void RefIntegral::transformZeroForm(const CellJacobianBatch& JVol,
  const Array<int>& isLocalFlag,  
  const double& coeff,
  RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(refIntegrationTimer());

  TEST_FOR_EXCEPTION(order() != 0, InternalError,
    "RefIntegral::transformZeroForm() called "
    "for form of order " << order());

  Tabs tabs;  
  SUNDANCE_MSG1(integrationVerb(), tabs << "doing zero form by reference");

  double& a = (*A)[0];
  int flops = 0;

  /* if we don't need to check whether elements are local, we
   * can streamline the loop. This will be the case when
   * we are evaluating a functional but not its gradient */
  double w = coeff * W_[0][0];
  if ((int) isLocalFlag.size()==0)
  {
    for (int c=0; c<JVol.numCells(); c++)
    {
      a += w * fabs(JVol.detJ()[c]);
    }
    flops = 2*JVol.numCells();
  }
  else
  {
    TEST_FOR_EXCEPTION( (int) isLocalFlag.size() != JVol.numCells(),
      RuntimeError,
      "mismatch between isLocalFlag.size()=" 
      << isLocalFlag.size()
      << " and JVol.numCells()=" << JVol.numCells());

    for (int c=0; c<JVol.numCells(); c++)
    {
      if (isLocalFlag[c]) 
      {
        flops+=2; 
        a += w * fabs(JVol.detJ()[c]);
      }
    }
  }
  addFlops(flops);
}

void RefIntegral::transformOneForm(const CellJacobianBatch& JTrans,  
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex,
  const double& coeff,
  RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(refIntegrationTimer());
  TEST_FOR_EXCEPTION(order() != 1, InternalError,
    "RefIntegral::transformOneForm() called for form "
    "of order " << order());
  Tabs tabs;  
  SUNDANCE_MSG1(integrationVerb(),
    tabs << "doing one form by reference");
  
  /* If the derivative order is zero, the only transformation to be done 
   * is to multiply by the cell's Jacobian determinant */
  if (testDerivOrder() == 0)
  {
    double* aPtr = &((*A)[0]);
    int count = 0;
    int nfc = nFacetCases();
    for (int c=0; c<JVol.numCells(); c++)
    {
      double detJ = coeff * fabs(JVol.detJ()[c]);
      int fc = 0;
      if (nfc != 1) fc = facetIndex[c];
      const Array<double>& w = W_[fc];
      for (int n=0; n<nNodes(); n++, count++) 
      {
        aPtr[count] += detJ*w[n];
      }
    }
    addFlops(JVol.numCells() * (nNodes() + 1));
  }
  else
  {
    /* If the derivative order is nonzero, then we have to do a transformation. 
     * If we're also on a cell of dimension lower than maximal, we need to refer
     * to the facet index of the facet being integrated. */
    int nCells = JVol.numCells();
    double one = 1.0;
    int nTransRows = nRefDerivTest();

    createOneFormTransformationMatrix(JTrans, JVol);

    SUNDANCE_MSG3(transformVerb(),
      Tabs() << "transformation matrix=" << G(alpha()));
    int nNodes0 = nNodes();
      
    if (nFacetCases()==1)
    {
      /* if we're on a maximal cell, we can do transformations 
       * for all cells in a single batch. 
       */
      ::dgemm_("N", "N", &nNodes0, &nCells, &nTransRows, &coeff, &(W_[0][0]),
        &nNodes0, &(G(alpha())[0]), &nTransRows, &one, 
        &((*A)[0]), &nNodes0);
    }
    else
    {
      /* If we're on a lower-dimensional cell and have to transform, 
       * we've got to do each transformation using a different facet case */
      int N = 1;
      for (int c=0; c<JVol.numCells(); c++)
      {
        int fc = facetIndex[c];
        double* aPtr = &((*A)[c*nNodes0]);
        ::dgemm_("N", "N", &nNodes0, &N, &nTransRows, &coeff, &(W_[fc][0]),
          &nNodes0, &(G(alpha())[c*nTransRows]), &nTransRows, &one, 
          aPtr, &nNodes0);
              
      }
          
    }
      
    addFlops(2 * nNodes0 * nCells * nTransRows);
  }
}

void RefIntegral::transformTwoForm(const CellJacobianBatch& JTrans,
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex, 
  const double& coeff,
  RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(refIntegrationTimer());
  TEST_FOR_EXCEPTION(order() != 2, InternalError,
    "RefIntegral::transformTwoForm() called for form "
    "of order " << order());
  
  Tabs tabs;  
  SUNDANCE_MSG1(transformVerb(), tabs << "doing two form by reference");

  /* If the derivative orders are zero, the only transformation to be done 
   * is to multiply by the cell's Jacobian determinant */
  if (testDerivOrder() == 0 && unkDerivOrder() == 0)
  {
    double* aPtr = &((*A)[0]);
    int count = 0;
    for (int c=0; c<JVol.numCells(); c++)
    {
      int fc = 0;
      if (nFacetCases() != 1) fc = facetIndex[c];
      const Array<double>& w = W_[fc];
      double detJ = coeff * fabs(JVol.detJ()[c]);
      for (int n=0; n<nNodes(); n++, count++) 
      {
        aPtr[count] += detJ*w[n];
      }
    }
    addFlops(JVol.numCells() * (nNodes() + 1));
  }
  else
  {
    /* If the derivative order is nonzero, then we have to do a transformation. 
     * If we're also on a cell of dimension lower than maximal, we need to refer
     * to the facet index of the facet being integrated. */
    int nCells = JVol.numCells();
    double one = 1.0;
    int nTransRows = nRefDerivUnk()*nRefDerivTest();

    createTwoFormTransformationMatrix(JTrans, JVol);
      
    double* GPtr;
    if (testDerivOrder() == 0)
    {
      GPtr = &(G(beta())[0]);
      SUNDANCE_MSG2(transformVerb(),
        Tabs() << "transformation matrix=" << G(beta()));
    }
    else if (unkDerivOrder() == 0)
    {
      GPtr = &(G(alpha())[0]);
      SUNDANCE_MSG2(transformVerb(),
        Tabs() << "transformation matrix=" << G(alpha()));
    }
    else
    {
      GPtr = &(G(alpha(), beta())[0]);
      SUNDANCE_MSG2(transformVerb(),
        Tabs() << "transformation matrix=" 
        << G(alpha(),beta()));
    }
      
    int nNodes0 = nNodes();

    if (nFacetCases()==1)
    {
      /* if we're on a maximal cell, we can do transformations 
       * for all cells in a single batch. 
       */
      ::dgemm_("N", "N", &nNodes0, &nCells, &nTransRows, &coeff, &(W_[0][0]),
        &nNodes0, GPtr, &nTransRows, &one, 
        &((*A)[0]), &nNodes0);
    }
    else
    {
      /* If we're on a lower-dimensional cell and have to transform, 
       * we've got to do each transformation using a different facet case */
      int N = 1;
      for (int c=0; c<JVol.numCells(); c++)
      {
        int fc = facetIndex[c];
        double* aPtr = &((*A)[c*nNodes0]);
        double* gPtr = &(GPtr[c*nTransRows]);
        SUNDANCE_MSG2(integrationVerb(),
          tabs << "c=" << c << ", facet case=" << fc
          << " W=" << W_[fc]);
        ::dgemm_("N", "N", &nNodes0, &N, &nTransRows, &coeff, &(W_[fc][0]),
          &nNodes0, gPtr, &nTransRows, &one, 
          aPtr, &nNodes0);
              
      }
          
    }
      
    addFlops(2 * nNodes0 * nCells * nTransRows);
  }
}
