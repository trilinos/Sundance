/* @HEADER@ */
/* @HEADER@ */

#include "SundanceRefBasisProductIntegral.hpp"
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



RefBasisProductIntegral:: RefBasisProductIntegral(int dim, 
                                                  const CellType& cellType,
                                                  const BasisFamily& testBasis,
                                                  int testOrder,
                                                  const QuadratureFamily& quad)
  : W_(),
    dim_(dim),
    testOrder_(testOrder), 
    nRefDerivTest_(ipow(dim, testOrder)),
    nNodesTest_(testBasis.nNodes(cellType)),
    unkOrder_(-1), 
    nRefDerivUnk_(-1),
    nNodesUnk_(-1),
    isTwoForm_(false)
{
  W_.resize(nRefDerivTest_ * nNodesTest_);
  for (int i=0; i<W_.size(); i++) W_[i]=0.0;

  Array<Array<Array<double> > > testBasisVals(nRefDerivTest());
        
  for (int r=0; r<nRefDerivTest(); r++)
    {
      MultiIndex mi;
      if (miTest.order()==1) mi[r] = 1;
      testBasis.ptr()->refEval(cellType, quadPts, mi, 
                               testBasisVals[r]);
    }
        
}



RefBasisProductIntegral::RefBasisProductIntegral(int dim,
                                                 const CellType& cellType,
                                                 const BasisFamily& testBasis,
                                                 int testOrder,
                                                 const BasisFamily& unkBasis,
                                                 int unkOrder,
                                                 const QuadratureFamily& quad)
  : W_(),
    dim_(dim),
    testOrder_(testOrder), 
    nRefDerivTest_(ipow(dim, testOrder)),
    nNodesTest_(testBasis.nNodes(cellType)), 
    unkOrder_(unkOrder), 
    nRefDerivUnk_(ipow(dim, unkOrder)),
    nNodesUnk_(unkBasis.nNodes(cellType)), 
    isTwoForm_(true)
{
  W_.resize(nRefDerivTest_ * nNodesTest_  * nRefDerivUnk_ * nNodesUnk_);
  for (int i=0; i<W_.size(); i++) W_[i]=0.0;

  Array<Array<Array<double> > > testBasisVals(nRefDerivTest());
  Array<Array<Array<double> > > unkBasisVals(nRefDerivUnk());
        
  for (int r=0; r<nRefDerivTest(); r++)
    {
      MultiIndex mi;
      if (miTest.order()==1) mi[r] = 1;
      testBasis.ptr()->refEval(cellType(), quadPts, mi, 
                               testBasisVals[r]);
    }
  for (int r=0; r<nRefDerivUnk(); r++)
    {
      MultiIndex mi;
      if (miUnk.order()==1) mi[r] = 1;
      unkBasis.ptr()->refEval(cellType(), quadPts, mi, unkBasisVals[r]);
    }

  for (int q=0; q<nQuad; q++)
    {
      for (int t=0; t<nRefDerivTest(); t++)
        {
          for (int nt=0; nt<nNodesTest(); nt++)
            {
              for (int u=0; u<nRefDerivTest(); u++)
                {
                  for (int nu=0; nu<nNodesUnk(); nu++)
                    {
                      value(t, nt, u, nu) 
                        += quadWgts[q] * testBasisVals[t][q][nt] 
                        * unkBasisVals[u][q][nu];
                    }
                }
            }
        }
    }
}

void RefBasisProductIntegral::transformTwoForm(const CellJacobianBatch& J,  
                                               const Array<int>& alpha,
                                               const Array<int>& beta,
                                               const Array<double>& coeff,
                                               RefCountPtr<Array<double> >& A) const
{
  /* If the derivative orders are zero, the only transformation to be done 
   * is to multiply by the cell's Jacobian determinant */
  if (testDerivOrder_ == 0 && unkDerivOrder_ == 0)
    {
      double* aPtr = &((*A)[0]);
      int count = 0;
      for (int c=0; c<J.numCells(); c++)
        {
          double detJ = coeff[0] * J.detJ()[c];
          for (int n=0; n<nNodes; n++, count++) 
            {
              aPtr[count] = detJ*W_[count];
            }
        }
    }
  else
    {
      createTwoFormTransformationMatrix(J, alpha, beta, coeff);
      
      ::dgemm_("N", "N", &nNodes, &nCells, &numTransRows, &one, &(W_[0]),
               &nNodes, &(G()[0]), &numTransRows, &zero, &((*A)[0]), &nCells);
       
    }
}

void RefBasisProductIntegral
::createTwoFormTransformationMatrix(const CellJacobianBatch& J,  
                                    const Array<int>& alpha,
                                    const Array<int>& beta,
                                    const Array<double>& coeff) const 
{
  

  /* If both derivative orders are 1, then we have to transform both
   * basis functions */
  
  if (testDerivOrder_ == 1 && unkDerivOrder_ == 1)
    {
      G().resize(J.numCells() * J.dim() * J.dim());

      for (int c=0; c<J.numCells(); j++)
        {
          double detJ = J.detJ()[c];
          for (int gamma=0; gamma<dim_; gamma++)
            {
              for (int delta=0; delta<dim_; delta++)
                {
                  double sum = 0.0;
                  for (int t=0; t<alpha.size(); t++)
                    {
                      sum += J.invJ(c, alpha[t], gamma)*J.invJ(c, beta[t], delta)
                        * coeff[t];
                    }
                  G()[c*dim_*dim_ + dim_*gamma + delta ] = detJ*sum; 
                }
            }
        }
    }

  else if (testDerivOrder_ == 1 && unkDerivOrder_ == 0)
    {
      G().resize(J.numCells() * J.dim());

      for (int c=0; c<J.numCells(); j++)
        {
          double detJ = J.detJ()[c];
          for (int gamma=0; gamma<dim; gamma++)
            {
              double sum = 0.0;
              for (int t=0; t<alpha.size(); t++)
                {
                  sum += J.invJ(c, alpha[t], gamma) * coeff[t];
                }
              G()[c*dim_ + gamma] = detJ*sum; 
            }
        }
    }
  else if (testDerivOrder_ == 0 && unkDerivOrder_ == 1)
    {
      G().resize(J.numCells() * J.dim());

      for (int c=0; c<J.numCells(); j++)
        {
          double detJ = J.detJ()[c];
          for (int delta=0; delta<dim; delta++)
            {
              double sum = 0.0;
              for (int t=0; t<beta.size(); t++)
                {
                  sum += J.invJ(c, beta[t], delta) * coeff[t];
                }
              G()[c*dim_ + delta] = detJ*sum; 
            }
        }
    }
  else
    {
      
    }
}
