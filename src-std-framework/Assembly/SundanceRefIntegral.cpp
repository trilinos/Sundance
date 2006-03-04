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
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

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

static Time& refIntegrationTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("ref integration"); 
  return *rtn;
}



RefIntegral::RefIntegral(int spatialDim,
                         const CellType& maxCellType,
                         int dim, 
                         const CellType& cellType)
  : ElementIntegral(spatialDim, maxCellType, dim, cellType), W_()
{
  Tabs tab0;
  verbosity() = classVerbosity();
  SUNDANCE_OUT(this->verbosity() > VerbSilent, 
               tab0 << "************* computing reference 0-form integrals ********" 
               << endl << tab0 << "cell type=" << cellType);
}

RefIntegral::RefIntegral(int spatialDim,
                         const CellType& maxCellType,
                         int dim, 
                         const CellType& cellType,
                         const BasisFamily& testBasis,
                         int alpha,
                         int testDerivOrder)
  : ElementIntegral(spatialDim, maxCellType, dim, cellType, 
                    testBasis, alpha, testDerivOrder), W_()
{
  Tabs tab0;
  verbosity() = classVerbosity();
  SUNDANCE_OUT(this->verbosity() > VerbSilent, 
               tab0 << "************* computing reference 1-form integrals ********" 
               << endl << tab0 << "test basis=" 
               << testBasis 
               << endl << tab0 << "cell type=" << cellType
               << endl << tab0 <<"differentiation order="
               << testDerivOrder);
  SUNDANCE_OUT(this->verbosity() > VerbMedium, 
               tab0 << "num test derivative components=" 
               << nRefDerivTest());

  TEST_FOR_EXCEPTION(testDerivOrder < 0 || testDerivOrder > 1,
                     InternalError,
                     "Test function derivative order=" << testDerivOrder
                     << " must be 0 or 1");

  W_.resize(nFacetCases());

  CellType evalCellType = cellType;
  if (nFacetCases() > 1) evalCellType = maxCellType;
  
  for (int fc=0; fc<nFacetCases(); fc++)
    {
      W_[fc].resize(nRefDerivTest() * nNodesTest());
      for (unsigned int i=0; i<W_[fc].size(); i++) W_[fc][i]=0.0;

      Array<Array<Array<double> > > testBasisVals(nRefDerivTest());
  
      QuadratureFamily quad = new GaussianQuadrature(max(1, testBasis.order()));
      Array<Point> quadPts;
      Array<double> quadWeights;
      if (nFacetCases()==1) quad.getPoints(cellType, quadPts, quadWeights);
      else quad.getFacetPoints(maxCellType, dim, fc, quadPts, quadWeights);

      int nQuad = quadPts.size();

      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "quad pts" << quadPts);

      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "quad weights" << quadWeights);

      for (int r=0; r<nRefDerivTest(); r++)
        {
          MultiIndex mi;
          if (testDerivOrder==1) mi[r] = 1;
          testBasis.ptr()->refEval(spatialDim, evalCellType, quadPts, mi, 
                                   testBasisVals[r]);
        }

      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "basis values" << testBasisVals);


      for (int q=0; q<nQuad; q++)
        {
          for (int t=0; t<nRefDerivTest(); t++)
            {
              for (int nt=0; nt<nNodesTest(); nt++)
                {
                  value(fc, t, nt) 
                    += chop(quadWeights[q] * testBasisVals[t][q][nt]) ;
                }
            }
        }    

      for (unsigned int i=0; i<W_[fc].size(); i++) W_[fc][i] = chop(W_[fc][i]);

      addFlops(3*nQuad*nRefDerivTest()*nNodesTest() + W_[fc].size());

    }
  if (verbosity() > VerbMedium)
    {
      print(cerr);
    }
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
                         int unkDerivOrder)
  : ElementIntegral(spatialDim, maxCellType,  dim, cellType,
                    testBasis, alpha, testDerivOrder, 
                    unkBasis, beta, unkDerivOrder), W_()

{
  Tabs tab0;
  verbosity() = classVerbosity();
  SUNDANCE_OUT(this->verbosity() > VerbSilent, 
               tab0 << " ************* computing reference 2-form integrals ***************" 
               << endl << tab0 << "test basis=" 
               << testBasis 
               << endl << tab0 << "unk basis=" << unkBasis
               << endl << tab0 << "cell type=" << cellType
               << endl << tab0 <<"differentiation order t=" 
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
      W_[fc].resize(nRefDerivTest() * nNodesTest()  * nRefDerivUnk() * nNodesUnk());
      for (unsigned int i=0; i<W_[fc].size(); i++) W_[fc][i]=0.0;

      Array<Array<Array<double> > > testBasisVals(nRefDerivTest());
      Array<Array<Array<double> > > unkBasisVals(nRefDerivUnk());
        
      QuadratureFamily quad 
        = new GaussianQuadrature(max(testBasis.order() + unkBasis.order(), 1));
      Array<Point> quadPts;
      Array<double> quadWeights;
      if (nFacetCases()==1) quad.getPoints(cellType, quadPts, quadWeights);
      else quad.getFacetPoints(maxCellType, dim, fc, quadPts, quadWeights);
      int nQuad = quadPts.size();

      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "quad pts" << quadPts);

      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "quad weights" << quadWeights);

      for (int r=0; r<nRefDerivTest(); r++)
        {
          MultiIndex mi;
          if (testDerivOrder==1) mi[r] = 1;
          testBasis.ptr()->refEval(spatialDim, evalCellType, quadPts, mi, 
                                   testBasisVals[r]);
        }

      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "test basis values" << testBasisVals);

      for (int r=0; r<nRefDerivUnk(); r++)
        {
          MultiIndex mi;
          if (unkDerivOrder==1) mi[r] = 1;
          unkBasis.ptr()->refEval(spatialDim, evalCellType, 
                                  quadPts, mi, unkBasisVals[r]);
        }


      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   tab0 << "unk basis values" << unkBasisVals);

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
                            += chop(quadWeights[q] * testBasisVals[t][q][nt] 
                                    * unkBasisVals[u][q][nu]);
                        }
                    }
                }
            }
        }
      addFlops(4*nQuad*nRefDerivTest()*nNodesTest()*nRefDerivUnk()*nNodesUnk()
               + W_[fc].size());
      for (unsigned int i=0; i<W_[fc].size(); i++) W_[fc][i] = chop(W_[fc][i]);
    }
  if (verbosity() > VerbMedium)
    {
      print(cerr);
    }

}


void RefIntegral::print(ostream& os) const 
{
  for (int fc=0; fc<nFacetCases(); fc++)
    {
      if (order()==1)
        {
          Tabs tab1;
          os << tab1 << "reference one-form values" << endl;
          if (testDerivOrder()==0)
            {
              Tabs tab2;
              os << tab2 << W_[fc] << endl;
            }
          else
            {
              Tabs tab2;
              for (int r=0; r<dim(); r++)
                {
                  os << tab2 << "dir=" << r << " {";
                  for (int n=0; n<nNodes(); n++) 
                    {
                      if (n != 0) os << ", ";
                      os << value(fc, r, n);
                    }
                  os << "}" << endl;
                }
            }
          os << tab1 << endl << tab1 << endl;
        }
      else if (order()==2)
        {
          Tabs tab1;
          os << tab1 << "reference two-form values" << endl;
          if (testDerivOrder()==0 && unkDerivOrder()==0)
            {
              Tabs tab2;
              os << tab2 << "{";
              for (int nt=0; nt<nNodesTest(); nt++) 
                {
                  if (nt!=0) os << ", ";
                  os << "{";
                  for (int nu=0; nu<nNodesUnk(); nu++)
                    {
                      if (nu!=0) os << ", ";
                      os << value(fc, 0, nt, 0, nu);
                    }
                  os << "}";
                }
              os << "}" << endl;
            }
          else if (testDerivOrder()==1 && unkDerivOrder()==1)
            {
              Tabs tab2;
              for (int t=0; t<dim(); t++)
                {
                  for (int u=0; u<dim(); u++)
                    {
                      os << tab2 << "test dir=" << t 
                         << ", unk dir=" << u << endl;
                      Tabs tab3;
                      os << tab3 << "{";
                      for (int nt=0; nt<nNodesTest(); nt++) 
                        {
                          if (nt!=0) os << ", ";
                          os << "{";
                          for (int nu=0; nu<nNodesUnk(); nu++)
                            {
                              if (nu!=0) os << ", ";
                              os << value(fc, t, nt, u, nu);
                            }
                          os << "}";
                        }
                      os << "}" << endl;
                    }
                }
            }
          else if (testDerivOrder()==1 && unkDerivOrder()==0)
            {
              Tabs tab2;
              for (int t=0; t<dim(); t++)
                {
                  os << tab2 << "test dir=" << t << endl;
                  Tabs tab3;
                  os << tab3 << "{";
                  for (int nt=0; nt<nNodesTest(); nt++) 
                    {
                      if (nt!=0) os << ", ";
                      os << "{";
                      for (int nu=0; nu<nNodesUnk(); nu++)
                        {
                          if (nu!=0) os << ", ";
                          os << value(fc, t, nt, 0, nu);
                        }
                      os << "}";
                    }
                  os << "}" << endl;
                }
            }
          else /* if (testDerivOrder()==1 && unkDerivOrder()==0) */
            {
              Tabs tab2;
              for (int u=0; u<dim(); u++)
                {
                  os << tab2 << "unk dir=" << u << endl;
                  Tabs tab3;
                  os << tab3 << "{";
                  for (int nt=0; nt<nNodesTest(); nt++) 
                    {
                      if (nt!=0) os << ", ";
                      os << "{";
                      for (int nu=0; nu<nNodesUnk(); nu++)
                        {
                          if (nu!=0) os << ", ";
                          os << value(fc, 0, nt, u, nu);
                        }
                      os << "}";
                    }
                  os << "}" << endl;
                }
            }
          os << tab1 << endl << tab1 << endl;
          os << tab1 << endl << tab1 << endl;
        }
    }
}


void RefIntegral::transformZeroForm(const CellJacobianBatch& JVol,  
                                    const double& coeff,
                                    RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(refIntegrationTimer());
  TEST_FOR_EXCEPTION(order() != 0, InternalError,
                     "RefIntegral::transformZeroForm() called "
                     "for form of order " << order());
  

  /* The result for each cell is the cell's Jacobian determinant */
  double& a = (*A)[0];
  for (int c=0; c<JVol.numCells(); c++)
    {
      a += coeff * fabs(JVol.detJ()[c]);
    }
  addFlops(2*JVol.numCells());
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
  
  /* If the derivative order is zero, the only transformation to be done 
   * is to multiply by the cell's Jacobian determinant */
  if (testDerivOrder() == 0)
    {
      double* aPtr = &((*A)[0]);
      int count = 0;
      const Array<double>& w = W_[0];
      for (int c=0; c<JVol.numCells(); c++)
        {
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
      int nTransRows = nRefDerivTest();

      createOneFormTransformationMatrix(JTrans, JVol);

      SUNDANCE_OUT(this->verbosity() > VerbMedium, 
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
  
  /* If the derivative orders are zero, the only transformation to be done 
   * is to multiply by the cell's Jacobian determinant */
  if (testDerivOrder() == 0 && unkDerivOrder() == 0)
    {
      double* aPtr = &((*A)[0]);
      const Array<double>& w = W_[0];
      int count = 0;
      for (int c=0; c<JVol.numCells(); c++)
        {
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
              ::dgemm_("N", "N", &nNodes0, &N, &nTransRows, &coeff, &(W_[fc][0]),
                       &nNodes0, gPtr, &nTransRows, &one, 
                       aPtr, &nNodes0);
              
            }
          
        }
      
      addFlops(2 * nNodes0 * nCells * nTransRows);
    }
}
