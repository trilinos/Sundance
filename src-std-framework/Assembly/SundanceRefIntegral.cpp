/* @HEADER@ */
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

RefIntegral::RefIntegral(int dim, 
                         const CellType& cellType,
                         const BasisFamily& testBasis,
                         int testDerivOrder)
  : W_(),
    dim_(dim),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(dim, testDerivOrder)),
    nNodesTest_(testBasis.nNodes(cellType)),
    unkDerivOrder_(-1), 
    nRefDerivUnk_(-1),
    nNodesUnk_(-1),
    nNodes_(nNodesTest_),
    isTwoForm_(false)
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
               << nRefDerivTest_);

  TEST_FOR_EXCEPTION(testDerivOrder < 0 || testDerivOrder > 1,
                     InternalError,
                     "Test function derivative order=" << testDerivOrder
                     << " must be 0 or 1");
  
  W_.resize(nRefDerivTest_ * nNodesTest_);
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
      if (testDerivOrder_==1) mi[r] = 1;
      testBasis.ptr()->refEval(cellType, quadPts, mi, 
                               testBasisVals[r]);
    }

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "basis values" << testBasisVals);


  for (int q=0; q<nQuad; q++)
    {
      for (int t=0; t<nRefDerivTest_; t++)
        {
          for (int nt=0; nt<nNodesTest_; nt++)
            {
              value(t, nt) 
                += chop(quadWeights[q] * testBasisVals[t][q][nt]) ;
            }
        }
    }    

  for (int i=0; i<W_.size(); i++) W_[i] = chop(W_[i]);

  if (verbosity() > VerbLow)
    {
      Tabs tab1;
      cerr << tab1 << "reference one-form values" << endl;
      if (testDerivOrder_==0)
        {
          Tabs tab2;
          cerr << tab2 << W_ << endl;
        }
      else
        {
          Tabs tab2;
          for (int r=0; r<dim_; r++)
            {
              cerr << tab2 << "dir=" << r << " {";
              for (int n=0; n<nNodes_; n++) 
                {
                  if (n != 0) cerr << ", ";
                  cerr << value(r, n);
                }
              cerr << "}" << endl;
            }
        }
      cerr << tab1 << endl << tab1 << endl;
    }
}



RefIntegral::RefIntegral(int dim,
                         const CellType& cellType,
                         const BasisFamily& testBasis,
                         int testDerivOrder,
                         const BasisFamily& unkBasis,
                         int unkDerivOrder)
  : W_(),
    dim_(dim),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(dim, testDerivOrder)),
    nNodesTest_(testBasis.nNodes(cellType)), 
    unkDerivOrder_(unkDerivOrder), 
    nRefDerivUnk_(ipow(dim, unkDerivOrder)),
    nNodesUnk_(unkBasis.nNodes(cellType)), 
    nNodes_(nNodesTest_*nNodesUnk_),
    isTwoForm_(true)
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
               << nRefDerivTest_);

  TEST_FOR_EXCEPTION(testDerivOrder < 0 || testDerivOrder > 1,
                     InternalError,
                     "Test function derivative order=" << testDerivOrder
                     << " must be 0 or 1");
  
  TEST_FOR_EXCEPTION(unkDerivOrder < 0 || unkDerivOrder > 1,
                     InternalError,
                     "Unknown function derivative order=" << unkDerivOrder
                     << " must be 0 or 1");

  
  W_.resize(nRefDerivTest_ * nNodesTest_  * nRefDerivUnk_ * nNodesUnk_);
  for (int i=0; i<W_.size(); i++) W_[i]=0.0;

  Array<Array<Array<double> > > testBasisVals(nRefDerivTest_);
  Array<Array<Array<double> > > unkBasisVals(nRefDerivUnk_);
        
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
      if (testDerivOrder_==1) mi[r] = 1;
      testBasis.ptr()->refEval(cellType, quadPts, mi, 
                               testBasisVals[r]);
    }

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "test basis values" << testBasisVals);

  for (int r=0; r<nRefDerivUnk(); r++)
    {
      MultiIndex mi;
      if (unkDerivOrder_) mi[r] = 1;
      unkBasis.ptr()->refEval(cellType, quadPts, mi, unkBasisVals[r]);
    }


  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "unk basis values" << unkBasisVals);

  for (int q=0; q<nQuad; q++)
    {
      for (int t=0; t<nRefDerivTest_; t++)
        {
          for (int nt=0; nt<nNodesTest_; nt++)
            {
              for (int u=0; u<nRefDerivUnk_; u++)
                {
                  for (int nu=0; nu<nNodesUnk_; nu++)
                    {
                      value(t, nt, u, nu) 
                        += chop(quadWeights[q] * testBasisVals[t][q][nt] 
                        * unkBasisVals[u][q][nu]);
                    }
                }
            }
        }
    }

  for (int i=0; i<W_.size(); i++) W_[i] = chop(W_[i]);

  if (verbosity() > VerbLow)
    {
      Tabs tab1;
      cerr << tab1 << "reference two-form values" << endl;
      if (testDerivOrder_==0 && unkDerivOrder_==0)
        {
          Tabs tab2;
          cerr << tab2 << "{";
          for (int nt=0; nt<nNodesTest_; nt++) 
            {
              if (nt!=0) cerr << ", ";
              cerr << "{";
              for (int nu=0; nu<nNodesUnk_; nu++)
                {
                  if (nu!=0) cerr << ", ";
                  cerr << value(0, nt, 0, nu);
                }
              cerr << "}";
            }
          cerr << "}" << endl;
        }
      else if (testDerivOrder_==1 && unkDerivOrder_==1)
        {
          Tabs tab2;
          for (int t=0; t<dim_; t++)
            {
              for (int u=0; u<dim_; u++)
                {
                  cerr << tab2 << "test dir=" << t 
                       << ", unk dir=" << u << endl;
                  Tabs tab3;
                  cerr << tab3 << "{";
                  for (int nt=0; nt<nNodesTest_; nt++) 
                    {
                      if (nt!=0) cerr << ", ";
                      cerr << "{";
                      for (int nu=0; nu<nNodesUnk_; nu++)
                        {
                          if (nu!=0) cerr << ", ";
                          cerr << value(t, nt, u, nu);
                        }
                      cerr << "}";
                    }
                  cerr << "}" << endl;
                }
            }
        }
      else if (testDerivOrder_==1 && unkDerivOrder_==0)
        {
          Tabs tab2;
          for (int t=0; t<dim_; t++)
            {
              cerr << tab2 << "test dir=" << t << endl;
              Tabs tab3;
              cerr << tab3 << "{";
              for (int nt=0; nt<nNodesTest_; nt++) 
                {
                  if (nt!=0) cerr << ", ";
                  cerr << "{";
                  for (int nu=0; nu<nNodesUnk_; nu++)
                    {
                      if (nu!=0) cerr << ", ";
                      cerr << value(t, nt, 0, nu);
                    }
                  cerr << "}" << endl;
                }
            }
        }
      else /* if (testDerivOrder_==1 && unkDerivOrder_==0) */
        {
          Tabs tab2;
          for (int u=0; u<dim_; u++)
            {
              cerr << tab2 << "unk dir=" << u << endl;
              Tabs tab3;
              cerr << tab3 << "{";
              for (int nt=0; nt<nNodesTest_; nt++) 
                {
                  if (nt!=0) cerr << ", ";
                  cerr << "{";
                  for (int nu=0; nu<nNodesUnk_; nu++)
                    {
                      if (nu!=0) cerr << ", ";
                      cerr << value(0, nt, u, nu);
                    }
                  cerr << "}" << endl;
                }

            }
        }
      cerr << tab1 << endl << tab1 << endl;
      cerr << tab1 << endl << tab1 << endl;
    }

}

int RefIntegral::ipow(int base, int power) 
{
  int rtn = 1;
  for (int i=0; i<power; i++) rtn *= base;
  return rtn;
}

void RefIntegral::transformTwoForm(const CellJacobianBatch& J,  
                                   const Array<int>& alpha,
                                   const Array<int>& beta,
                                   const Array<double>& coeff,
                                   RefCountPtr<Array<double> >& A) const
{
  /* If the derivative orders are zero, the only transformation to be done 
   * is to multiply by the cell's Jacobian determinant */
  if (testDerivOrder_ == 0 && unkDerivOrder_ == 0)
    {
      A->resize(J.numCells() * nNodes_);
      double* aPtr = &((*A)[0]);
      int count = 0;
      for (int c=0; c<J.numCells(); c++)
        {
          double detJ = coeff[0] * J.detJ()[c];
          for (int n=0; n<nNodes_; n++, count++) 
            {
              aPtr[count] = detJ*W_[n];
            }
        }
    }
  else
    {
      int nCells = J.numCells();
      double one = 1.0;
      double zero = 0.0;
      int nTransRows = nRefDerivUnk_*nRefDerivTest_;
      A->resize(J.numCells() * nNodes_); 
      int info=0;

      createTwoFormTransformationMatrix(J, alpha, beta, coeff);
      
      ::dgemm_("N", "N", &nNodes_, &nCells, &nTransRows, &one, &(W_[0]),
               &nNodes_, &(G()[0]), &nTransRows, &zero, &((*A)[0]), &nNodes_);
       
    }
}

void RefIntegral
::createTwoFormTransformationMatrix(const CellJacobianBatch& J,  
                                    const Array<int>& alpha,
                                    const Array<int>& beta,
                                    const Array<double>& coeff) const 
{
  TEST_FOR_EXCEPTION(J.cellDim() != dim_, InternalError,
                     "Inconsistency between Jacobian dimension " << J.cellDim()
                     << " and cell dimension " << dim_ 
                     << " in RefIntegral::createTwoFormTransformationMatrix()");

  /* If both derivative orders are 1, then we have to transform both
   * basis functions */
  
  if (testDerivOrder_ == 1 && unkDerivOrder_ == 1)
    {
      G().resize(J.numCells() * J.cellDim() * J.cellDim());

      for (int c=0; c<J.numCells(); c++)
        {
          double detJ = J.detJ()[c];
          Array<double> invJ;
          J.getInvJ(c, invJ);
          for (int gamma=0; gamma<dim_; gamma++)
            {
              for (int delta=0; delta<dim_; delta++)
                {
                  double sum = 0.0;
                  for (int t=0; t<alpha.size(); t++)
                    {
                      sum += invJ[alpha[t] + gamma*dim_]
                        * invJ[beta[t]*dim_+delta]
                        * coeff[t];
                    }
                  G()[c*dim_*dim_ + dim_*gamma + delta ] = detJ*sum; 
                }
            }
        }
    }

  else if (testDerivOrder_ == 1 && unkDerivOrder_ == 0)
    {
      G().resize(J.numCells() * J.cellDim());

      for (int c=0; c<J.numCells(); c++)
        {
          double detJ = J.detJ()[c];
          Array<double> invJ;
          J.getInvJ(c, invJ);
          for (int gamma=0; gamma<dim_; gamma++)
            {
              double sum = 0.0;
              for (int t=0; t<alpha.size(); t++)
                {
                  sum += invJ[alpha[t] + dim_ * gamma] * coeff[t];
                }
              G()[c*dim_ + gamma] = detJ*sum; 
            }
        }
    }
  else /* if (testDerivOrder_ == 0 && unkDerivOrder_ == 1) */
    {
      G().resize(J.numCells() * J.cellDim());

      for (int c=0; c<J.numCells(); c++)
        {
          double detJ = J.detJ()[c];
          Array<double> invJ;
          J.getInvJ(c, invJ);
          for (int delta=0; delta<dim_; delta++)
            {
              double sum = 0.0;
              for (int t=0; t<beta.size(); t++)
                {
                  sum += invJ[beta[t] + dim_ * delta] * coeff[t];
                }
              G()[c*dim_ + delta] = detJ*sum; 
            }
        }
    }
}
