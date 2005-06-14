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



RefIntegral::RefIntegral(int dim, 
                         const CellType& cellType)
  : ElementIntegral(dim, cellType), W_()
{
  Tabs tab0;
  verbosity() = classVerbosity();
  SUNDANCE_OUT(verbosity() > VerbSilent, 
               tab0 << "************* computing reference 0-form integrals ********" 
               << endl << tab0 << "cell type=" << cellType);
}

RefIntegral::RefIntegral(int dim, 
                         const CellType& cellType,
                         const BasisFamily& testBasis,
                         int alpha,
                         int testDerivOrder)
  : ElementIntegral(dim, cellType, testBasis, alpha, testDerivOrder), W_()
{
  Tabs tab0;
  verbosity() = classVerbosity();
  SUNDANCE_OUT(verbosity() > VerbSilent, 
               tab0 << "************* computing reference 1-form integrals ********" 
               << endl << tab0 << "test basis=" 
               << testBasis 
               << endl << tab0 << "cell type=" << cellType
               << endl << tab0 <<"differentiation order="
               << testDerivOrder);
  SUNDANCE_OUT(verbosity() > VerbMedium, 
               tab0 << "num test derivative components=" 
               << nRefDerivTest());

  TEST_FOR_EXCEPTION(testDerivOrder < 0 || testDerivOrder > 1,
                     InternalError,
                     "Test function derivative order=" << testDerivOrder
                     << " must be 0 or 1");
  
  W_.resize(nRefDerivTest() * nNodesTest());
  for (int i=0; i<W_.size(); i++) W_[i]=0.0;

  Array<Array<Array<double> > > testBasisVals(nRefDerivTest());
  
  QuadratureFamily quad = new GaussianQuadrature(max(1, testBasis.order()));
  Array<Point> quadPts;
  Array<double> quadWeights;
  quad.getPoints(cellType, quadPts, quadWeights);
  int nQuad = quadPts.size();

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "quad pts" << quadPts);

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "quad weights" << quadWeights);

  for (int r=0; r<nRefDerivTest(); r++)
    {
      MultiIndex mi;
      if (testDerivOrder==1) mi[r] = 1;
      testBasis.ptr()->refEval(cellType, quadPts, mi, 
                               testBasisVals[r]);
    }

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "basis values" << testBasisVals);


  for (int q=0; q<nQuad; q++)
    {
      for (int t=0; t<nRefDerivTest(); t++)
        {
          for (int nt=0; nt<nNodesTest(); nt++)
            {
              value(t, nt) 
                += chop(quadWeights[q] * testBasisVals[t][q][nt]) ;
            }
        }
    }    

  for (int i=0; i<W_.size(); i++) W_[i] = chop(W_[i]);

  addFlops(3*nQuad*nRefDerivTest()*nNodesTest() + W_.size());

  if (verbosity() > VerbMedium)
    {
      print(cerr);
    }
}




RefIntegral::RefIntegral(int dim,
                         const CellType& cellType,
                         const BasisFamily& testBasis,
                         int alpha,
                         int testDerivOrder,
                         const BasisFamily& unkBasis,
                         int beta,
                         int unkDerivOrder)
  : ElementIntegral(dim, cellType, 
                    testBasis, alpha, testDerivOrder, 
                    unkBasis, beta, unkDerivOrder), W_()

{
  Tabs tab0;
  verbosity() = classVerbosity();
  SUNDANCE_OUT(verbosity() > VerbSilent, 
               tab0 << " ************* computing reference 2-form integrals ***************" 
               << endl << tab0 << "test basis=" 
               << testBasis 
               << endl << tab0 << "unk basis=" << unkBasis
               << endl << tab0 << "cell type=" << cellType
               << endl << tab0 <<"differentiation order t=" 
               << testDerivOrder << ", u=" << unkDerivOrder);
  SUNDANCE_OUT(verbosity() > VerbMedium, 
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

  
  W_.resize(nRefDerivTest() * nNodesTest()  * nRefDerivUnk() * nNodesUnk());
  for (int i=0; i<W_.size(); i++) W_[i]=0.0;

  Array<Array<Array<double> > > testBasisVals(nRefDerivTest());
  Array<Array<Array<double> > > unkBasisVals(nRefDerivUnk());
        
  QuadratureFamily quad 
    = new GaussianQuadrature(max(testBasis.order() + unkBasis.order(), 1));
  Array<Point> quadPts;
  Array<double> quadWeights;
  quad.getPoints(cellType, quadPts, quadWeights);
  int nQuad = quadPts.size();

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "quad pts" << quadPts);

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "quad weights" << quadWeights);

  for (int r=0; r<nRefDerivTest(); r++)
    {
      MultiIndex mi;
      if (testDerivOrder==1) mi[r] = 1;
      testBasis.ptr()->refEval(cellType, quadPts, mi, 
                               testBasisVals[r]);
    }

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "test basis values" << testBasisVals);

  for (int r=0; r<nRefDerivUnk(); r++)
    {
      MultiIndex mi;
      if (unkDerivOrder==1) mi[r] = 1;
      unkBasis.ptr()->refEval(cellType, quadPts, mi, unkBasisVals[r]);
    }


  SUNDANCE_OUT(verbosity() > VerbHigh, 
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
                      value(t, nt, u, nu) 
                        += chop(quadWeights[q] * testBasisVals[t][q][nt] 
                                * unkBasisVals[u][q][nu]);
                    }
                }
            }
        }
    }
  addFlops(4*nQuad*nRefDerivTest()*nNodesTest()*nRefDerivUnk()*nNodesUnk()
           + W_.size());
  for (int i=0; i<W_.size(); i++) W_[i] = chop(W_[i]);

  if (verbosity() > VerbMedium)
    {
      print(cerr);
    }

}


void RefIntegral::print(ostream& os) const 
{
  if (order()==1)
    {
      Tabs tab1;
      os << tab1 << "reference one-form values" << endl;
      if (testDerivOrder()==0)
        {
          Tabs tab2;
          os << tab2 << W_ << endl;
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
                  os << value(r, n);
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
                  os << value(0, nt, 0, nu);
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
                          os << value(t, nt, u, nu);
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
                      os << value(t, nt, 0, nu);
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
                      os << value(0, nt, u, nu);
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


void RefIntegral::transformZeroForm(const CellJacobianBatch& J,  
                                    const double& coeff,
                                    RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(refIntegrationTimer());
  TEST_FOR_EXCEPTION(order() != 0, InternalError,
                     "RefIntegral::transformZeroForm() called "
                     "for form of order " << order());
  

  /* The result for each cell is the cell's Jacobian determinant */
  double& a = (*A)[0];
  for (int c=0; c<J.numCells(); c++)
    {
      a += coeff * fabs(J.detJ()[c]);
    }
  addFlops(2*J.numCells());
}

void RefIntegral::transformOneForm(const CellJacobianBatch& J,  
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
      for (int c=0; c<J.numCells(); c++)
        {
          double detJ = coeff * fabs(J.detJ()[c]);
          for (int n=0; n<nNodes(); n++, count++) 
            {
              aPtr[count] += detJ*W_[n];
            }
        }
      addFlops(J.numCells() * (nNodes() + 1));
    }
  else
    {
      int nCells = J.numCells();
      double one = 1.0;
      double zero = 0.0;
      int nTransRows = nRefDerivTest();
      int info=0;

      createOneFormTransformationMatrix(J);

      SUNDANCE_OUT(verbosity() > VerbMedium, 
                   Tabs() << "transformation matrix=" << G(alpha()));
      
      int nNodes0 = nNodes();
      ::dgemm_("N", "N", &nNodes0, &nCells, &nTransRows, &coeff, &(W_[0]),
               &nNodes0, &(G(alpha())[0]), &nTransRows, &one, 
               &((*A)[0]), &nNodes0);

      addFlops(2 * nNodes0 * nCells * nTransRows);
    }
}

void RefIntegral::transformTwoForm(const CellJacobianBatch& J,  
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
      //      A->resize(J.numCells() * nNodes());
      double* aPtr = &((*A)[0]);
      int count = 0;
      for (int c=0; c<J.numCells(); c++)
        {
          double detJ = coeff * fabs(J.detJ()[c]);
          for (int n=0; n<nNodes(); n++, count++) 
            {
              aPtr[count] += detJ*W_[n];
            }
        }
      addFlops(J.numCells() * (nNodes() + 1));
    }
  else
    {
      int nCells = J.numCells();
      double one = 1.0;
      double zero = 0.0;
      int nTransRows = nRefDerivUnk()*nRefDerivTest();
      //      A->resize(J.numCells() * nNodes()); 
      int info=0;

      createTwoFormTransformationMatrix(J);
      
      double* GPtr;
      if (testDerivOrder() == 0)
        {
          GPtr = &(G(beta())[0]);
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       Tabs() << "transformation matrix=" << G(beta()));
        }
      else if (unkDerivOrder() == 0)
        {
          GPtr = &(G(alpha())[0]);
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       Tabs() << "transformation matrix=" << G(alpha()));
        }
      else
        {
          GPtr = &(G(alpha(), beta())[0]);
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       Tabs() << "transformation matrix=" 
                       << G(alpha(),beta()));
        }
      
      int nNodes0 = nNodes();
      ::dgemm_("N", "N", &nNodes0, &nCells, &nTransRows, &coeff, &(W_[0]),
               &nNodes0, GPtr, 
               &nTransRows, &one, &((*A)[0]), &nNodes0);
       
      addFlops(2 * nNodes0 * nCells * nTransRows);
    }
}

