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


void RefBasisProductIntegral::transformTwoForm(const CellJacobianBatch& J, 
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
          double detJ = J.detJ()[c];
          for (int n=0; n<nNodes; n++, count++) 
            {
              aPtr[count] = detJ*W_[count];
            }
        }
    }
  else
    {
      createTwoFormTransformationMatrix();
      
      ::dgemm_("N", "N", &nNodes, &nCells, &numTransRows, &one, &(W_[0]),
               &nNodes, &(G()[0]), &numTransRows, &zero, &((*A)[0]), &nCells);
       
    }
}

void RefBasisProductIntegral::createTwoFormTransformationMatrix() const 
{
  

  /* If both derivative orders are 1, then we have to transform both
   * basis functions */
  
  if (testDerivOrder_ == 1 && unkDerivOrder_ == 1)
    {
      G().resize(J.numCells() * J.dim() * J.dim());

      for (int c=0; c<J.numCells(); j++)
        {
          double detJ = J.detJ()[c];
          for (int gamma=0; gamma<dim; gamma++)
            {
              for (int delta=0; delta<dim; delta++)
                {
                  double sum = 0.0;
                  for (int t=0; t<alpha_.size(); t++)
                    {
                      sum += J.invJ(c, alpha_[t], gamma)*J.invJ(c, beta_[t], delta)
                            * coeff_[t];
                    }
                  G()[c*dim*dim + dim*gamma + delta ] = detJ*sum; 
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
              for (int t=0; t<alpha_.size(); t++)
                {
                  sum += J.invJ(c, alpha_[t], gamma) * coeff_[t];
                }
              G()[c*dim + gamma] = detJ*sum; 
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
              for (int t=0; t<beta_.size(); t++)
                {
                  sum += J.invJ(c, beta_[t], delta) * coeff_[t];
                }
              G()[c*dim + delta] = detJ*sum; 
            }
        }
    }
  else
    {
      
    }
}
