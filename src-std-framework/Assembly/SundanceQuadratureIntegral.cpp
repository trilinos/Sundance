/* @HEADER@ */
/* @HEADER@ */

#include "SundanceQuadratureIntegral.hpp"
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

static Time& quadratureTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("quadrature"); 
  return *rtn;
}


static Time& quadTransCreationTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("building transformation matrices for quad"); 
  return *rtn;
}


QuadratureIntegral::QuadratureIntegral(int dim, 
                                       const CellType& cellType,
                                       const QuadratureFamily& quad)
  : ElementIntegral(dim, cellType),
    W_(),
    nQuad_(0),
    useSumFirstMethod_(true)
{
  Tabs tab0;
  verbosity() = classVerbosity();

  /* create the quad points and weights */
  Array<double> quadWeights;
  Array<Point> quadPts;
  quad.getPoints(cellType, quadPts, quadWeights);
  nQuad_ = quadPts.size();
  
  W_.resize(nQuad());

  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab0 << "num quad pts" << nQuad());

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "quad weights" << quadWeights);

  for (int q=0; q<nQuad(); q++)
    {
      W_[q] = quadWeights[q];
    }    
}


QuadratureIntegral::QuadratureIntegral(int dim, 
                                       const CellType& cellType,
                                       const BasisFamily& testBasis,
                                       const Array<int>& alpha,
                                       int testDerivOrder,
                                       const QuadratureFamily& quad)
  : ElementIntegral(dim, cellType, testBasis, alpha, testDerivOrder),
    W_(),
    nQuad_(0),
    useSumFirstMethod_(true)
{
  Tabs tab0;
  verbosity() = classVerbosity();
  SUNDANCE_OUT(verbosity() > VerbSilent, 
               tab0 
               << "******** computing basis functions on quad pts *******"
               << endl << tab0 << "test basis=" 
               << testBasis 
               << endl << tab0 << "cell type=" << cellType
               << endl << tab0 << "differentiation order="
               << testDerivOrder);
  SUNDANCE_OUT(verbosity() > VerbMedium, 
               tab0 << "num test derivative components=" 
               << nRefDerivTest());

  TEST_FOR_EXCEPTION(testDerivOrder < 0 || testDerivOrder > 1,
                     InternalError,
                     "Test function derivative order=" << testDerivOrder
                     << " must be 0 or 1");
  

  /* create the quad points and weights */
  Array<double> quadWeights;
  Array<Point> quadPts;
  quad.getPoints(cellType, quadPts, quadWeights);
  nQuad_ = quadPts.size();
  
  W_.resize(nQuad() * nRefDerivTest() * nNodesTest());

  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab0 << "num quad pts" << nQuad());

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "quad pts" << quadPts);

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "quad weights" << quadWeights);

  Array<Array<Array<double> > > testBasisVals(nRefDerivTest());

  for (int r=0; r<nRefDerivTest(); r++)
    {
      MultiIndex mi;
      if (testDerivOrder==1) mi[r] = 1;
      testBasis.ptr()->refEval(cellType, quadPts, mi, 
                               testBasisVals[r]);
    }

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "basis values" << testBasisVals);


  for (int q=0; q<nQuad(); q++)
    {
      for (int t=0; t<nRefDerivTest(); t++)
        {
          for (int nt=0; nt<nNodesTest(); nt++)
            {
              wValue(q, t, nt) 
                = chop(quadWeights[q] * testBasisVals[t][q][nt]) ;
            }
        }
    }    
  
  addFlops(2*nQuad()*nRefDerivTest()*nNodesTest());

  if (verbosity() > VerbMedium)
    {
      print(cerr);
    }
}




QuadratureIntegral::QuadratureIntegral(int dim,
                                       const CellType& cellType,
                                       const BasisFamily& testBasis,
                                       const Array<int>& alpha,
                                       int testDerivOrder,
                                       const BasisFamily& unkBasis,
                                       const Array<int>& beta,
                                       int unkDerivOrder,
                                       const QuadratureFamily& quad)
  : ElementIntegral(dim, cellType, 
                    testBasis, alpha, testDerivOrder, 
                    unkBasis, beta, unkDerivOrder), 
    W_(),
    nQuad_(0),
    useSumFirstMethod_(true)
{
  Tabs tab0;
  verbosity() = classVerbosity();
  SUNDANCE_OUT(verbosity() > VerbSilent, 
               tab0 << " ************* computing basis func products on quad pts ***************" 
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

  /* get the quad pts and weights */
  Array<double> quadWeights;
  Array<Point> quadPts;
  quad.getPoints(cellType, quadPts, quadWeights);
  nQuad_ = quadPts.size();

  W_.resize(nQuad() * nRefDerivTest() * nNodesTest()  
            * nRefDerivUnk() * nNodesUnk());

  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab0 << "num quad pts" << nQuad());

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "quad pts" << quadPts);

  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab0 << "quad weights" << quadWeights);



  /* compute the basis functions */
  Array<Array<Array<double> > > testBasisVals(nRefDerivTest());
  Array<Array<Array<double> > > unkBasisVals(nRefDerivUnk());

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
                      wValue(q, t, nt, u, nu)
                        = chop(quadWeights[q] * testBasisVals[t][q][nt] 
                               * unkBasisVals[u][q][nu]);
                    }
                }
            }
        }
    }

  addFlops(3*nQuad()*nRefDerivTest()*nNodesTest()*nRefDerivUnk()*nNodesUnk()
           + W_.size());
  for (int i=0; i<W_.size(); i++) W_[i] = chop(W_[i]);

  if (verbosity() > VerbMedium)
    {
      print(cerr);
    }

}


void QuadratureIntegral::print(ostream& os) const 
{
  
}

void QuadratureIntegral::transformZeroForm(const CellJacobianBatch& J,  
                                           const double* const coeff,
                                           RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(quadratureTimer());
  TEST_FOR_EXCEPTION(order() != 0, InternalError,
                     "QuadratureIntegral::transformZeroForm() called "
                     "for form of order " << order());

  int flops = 0;

  A->resize(1);
  double& a = (*A)[0];
  a = 0.0;
  double* coeffPtr = (double*) coeff;
  for (int c=0; c<J.numCells(); c++)
    {
      double detJ = fabs(J.detJ()[c]);
      for (int q=0; q<nQuad(); q++, coeffPtr++)
        {
          a += W_[q]*(*coeffPtr)*detJ;
        }
    }

  addFlops(J.numCells()*(1 + 2*nQuad()));
}



void QuadratureIntegral::transformOneForm(const CellJacobianBatch& J,  
                                          const double* const coeff,
                                          RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(quadratureTimer());
  TEST_FOR_EXCEPTION(order() != 1, InternalError,
                     "QuadratureIntegral::transformOneForm() called for form "
                     "of order " << order());

  int flops = 0;

  /* If the derivative order is zero, the only thing to be done 
   * is to multiply by the cell's Jacobian determinant and sum over the
   * quad points */
  if (testDerivOrder() == 0)
    {
      A->resize(J.numCells() * nNodes());
      double* aPtr = &((*A)[0]);
      for (int i=0; i<A->size(); i++) aPtr[i] = 0.0;
      double* coeffPtr = (double*) coeff;
      int offset = 0 ;

      for (int c=0; c<J.numCells(); c++, offset+=nNodes())
        {
          double detJ = fabs(J.detJ()[c]);
          for (int q=0; q<nQuad(); q++, coeffPtr++)
            {
              double f = (*coeffPtr)*detJ;
              for (int n=0; n<nNodes(); n++) 
                {
                  aPtr[offset+n] += f*W_[n + nNodes()*q];
                }
            }
        }
      addFlops( J.numCells() * (1 + nQuad() * (1 + 2*nNodes())) );
    }
  else
    {
      int nCells = J.numCells();
      A->resize(nCells * nNodes()); 

      createOneFormTransformationMatrix(J, alpha());

      SUNDANCE_OUT(verbosity() > VerbMedium, 
                   Tabs() << "transformation matrix=" << G());
      
      if (useSumFirstMethod())
        {
          transformSummingFirst(J.numCells(), coeff, A);
        }
      else
        {
          transformSummingLast(J.numCells(), coeff, A);
        }
    }
  addFlops(flops);
}


void QuadratureIntegral::transformTwoForm(const CellJacobianBatch& J,  
                                          const double* const coeff,
                                          RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(quadratureTimer());
  TEST_FOR_EXCEPTION(order() != 2, InternalError,
                     "QuadratureIntegral::transformTwoForm() called for form "
                     "of order " << order());

  int flops = 0 ;

  /* If the derivative orders are zero, the only thing to be done 
   * is to multiply by the cell's Jacobian determinant and sum over the
   * quad points */
  if (testDerivOrder() == 0 && unkDerivOrder() == 0)
    {

      A->resize(J.numCells() * nNodes());
      double* aPtr = &((*A)[0]);
      for (int i=0; i<A->size(); i++) aPtr[i] = 0.0;
      double* coeffPtr = (double*) coeff;
      int offset = 0 ;

      for (int c=0; c<J.numCells(); c++, offset+=nNodes())
        {
          double detJ = fabs(J.detJ()[c]);
          for (int q=0; q<nQuad(); q++, coeffPtr++)
            {
              double f = (*coeffPtr)*detJ;
              for (int n=0; n<nNodes(); n++) 
                {
                  aPtr[offset+n] += f*W_[n + nNodes()*q];
                }
            }
        }

      addFlops( J.numCells() * (1 + nQuad() * (1 + 2*nNodes())) );
    }
  else
    {
      int nCells = J.numCells();
      A->resize(nCells * nNodes()); 

      createTwoFormTransformationMatrix(J, alpha(), beta());

      SUNDANCE_OUT(verbosity() > VerbMedium, 
                   Tabs() << "transformation matrix=" << G());
      
      if (useSumFirstMethod())
       {
         transformSummingFirst(J.numCells(), coeff, A);
       }
      else
        {
          transformSummingLast(J.numCells(), coeff, A);
        }
    }
  addFlops(flops);
}

void QuadratureIntegral
::transformSummingFirst(int nCells,
                        const double* const coeff,
                        RefCountPtr<Array<double> >& A) const
{
  double* aPtr = &((*A)[0]);
  double* coeffPtr = (double*) coeff;
  
  int nAdds = 0;
  int nMults = 0;

  //  cerr << "transform summing first " << endl;

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
      cerr << "nCells = " << nCells << endl;
      cerr << "nNodes = " << nNodes() << endl;
      cerr << "nQuad = " << nQuad() << endl;
      cerr << "transSize = " << transSize << endl;
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

      for (int q=0; q<nQuad(); q++, coeffPtr++)
        {
          double f = (*coeffPtr);
          //          cerr << "cell=" << c << " q=" << q << " f=" << f << endl;
          for (int n=0; n<swSize; n++, nMults++, nAdds++) 
            {
              sumWorkspace[n] += f*W_[n + q*swSize];
            }
        }
      //      cerr << "cell=" << c << " untransformed sum=" << sumWorkspace << endl;
      /* transform the sum */
      double* gCell = &(G()[transSize*c]);
      double* aCell = aPtr + nNodes()*c;
      for (int i=0; i<nNodes(); i++)
        {
          aCell[i] = 0.0;
          for (int j=0; j<transSize; j++, nMults++, nAdds++)
            {
              aCell[i] += sumWorkspace[nNodes()*j + i] * gCell[j];
            }
          //          cerr << "cell=" << c << " transformed sum= " << aCell[i] << endl;
        }
    }
  
  int flops = 2*(nCells * nNodes() * transSize) * (nQuad() + 1) ;
//   cerr << "number of multiplies: " << nMults << endl;
//   cerr << "number of adds: " << nAdds << endl;
//   cerr << "predicted flops = " << flops << endl;
//   cerr << "counted flops = " << nMults + nAdds << endl;
  addFlops(flops);
}

void QuadratureIntegral
::transformSummingLast(int nCells,
                       const double* const coeff,
                       RefCountPtr<Array<double> >& A) const
{
  double* aPtr = &((*A)[0]);
  int transSize = 0; 
  
  int nAdds = 0;
  int nMults = 0;

  //  cerr << "transform summing last " << endl;

  if (order()==2)
    {
      transSize = nRefDerivTest() * nRefDerivUnk();
    }
  else
    {
      transSize = nRefDerivTest();
    }


//   if (true)
//     {
//       cerr << "nCells = " << nCells << endl;
//       cerr << "nNodes = " << nNodes() << endl;
//       cerr << "nQuad = " << nQuad() << endl;
//       cerr << "transSize = " << transSize << endl;
//     }

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
      double* gCell = &(G()[transSize*c]);
      double* aCell = aPtr + nNodes()*c;
      for (int n=0; n<nNodes(); n++) aCell[n] = 0.0;

      for (int q=0; q<nQuad(); q++)
        {
          double f = coeff[c*nQuad() + q];
          //          cerr << "cell=" << c << " q=" << q << " f=" << f << endl;
          for (int t=0; t<transSize; t++, nMults++) jWorkspace[t]=f*gCell[t];

          for (int n=0; n<nNodes(); n++)
            {
              for (int t=0; t<transSize; t++, nMults++, nAdds++)
                {
                  aCell[n] += jWorkspace[t]*W_[n + nNodes()*(t + transSize*q)];
                }
            }
        }
    }
  int flops = nCells * nQuad() * transSize * (1 + 2*nNodes()) ;
//   cerr << "number of multiplies: " << nMults << endl;
//   cerr << "number of adds: " << nAdds << endl;
//   cerr << "predicted flops = " << flops << endl;
//   cerr << "counted flops = " << nMults + nAdds << endl;
  addFlops(flops);
}

void QuadratureIntegral
::createTwoFormTransformationMatrix(const CellJacobianBatch& J,  
                                    const Array<int>& alpha,
                                    const Array<int>& beta) const 
{
  TimeMonitor timer(quadTransCreationTimer());

  TEST_FOR_EXCEPTION(J.cellDim() != dim(), InternalError,
                     "Inconsistency between Jacobian dimension " << J.cellDim()
                     << " and cell dimension " << dim() 
                     << " in QuadratureIntegral::createTwoFormTransformationMatrix()");

  /* If both derivative orders are 1, then we have to transform both
   * basis functions */

  int flops = 0;
  
  if (testDerivOrder() == 1 && unkDerivOrder() == 1)
    {
      G().resize(J.numCells() * J.cellDim() * J.cellDim());

      flops = J.numCells() * dim() * (1 + 2*alpha.size());

      SUNDANCE_OUT(verbosity() > VerbMedium, 
                   Tabs() << "both derivs are first order");
      for (int c=0; c<J.numCells(); c++)
        {
          static Array<double> invJ;
          J.getInvJ(c, invJ);
          //          cerr << "cell=" << c << " invJ=" << invJ << endl;
          double detJ = fabs(J.detJ()[c]);
          for (int gamma=0; gamma<dim(); gamma++)
            {
              for (int delta=0; delta<dim(); delta++)
                {
                  double sum = 0.0;
                  for (int t=0; t<alpha.size(); t++)
                    {
                      sum += invJ[alpha[t] + gamma*dim()]
                        * invJ[beta[t]+ dim()*delta];
                    }
                  G()[dim()*(c*dim() + gamma) + delta ] = detJ*sum; 
                }
            }
        }
      //      cerr << "G = " << G() << endl;
    }

  else if (testDerivOrder() == 1 && unkDerivOrder() == 0)
    {
      G().resize(J.numCells() * J.cellDim());

      flops = J.numCells() * dim() * (1 + alpha.size());
      for (int c=0; c<J.numCells(); c++)
        {
          static Array<double> invJ;
          J.getInvJ(c, invJ);
          double detJ = fabs(J.detJ()[c]);
          for (int gamma=0; gamma<dim(); gamma++)
            {
              double sum = 0.0;
              for (int t=0; t<alpha.size(); t++)
                {
                  sum += invJ[alpha[t] + dim() * gamma];
                }
              G()[c*dim() + gamma] = detJ*sum; 
            }
        }
    }
  else /* if (testDerivOrder() == 0 && unkDerivOrder() == 1) */
    {
      G().resize(J.numCells() * J.cellDim());

      flops = J.numCells() * dim() * (1 + alpha.size());

      for (int c=0; c<J.numCells(); c++)
        {
          static Array<double> invJ;
          J.getInvJ(c, invJ);
          double detJ = fabs(J.detJ()[c]);
          for (int delta=0; delta<dim(); delta++)
            {
              double sum = 0.0;
              for (int t=0; t<beta.size(); t++)
                {
                  sum += invJ[beta[t] + dim() * delta];
                }
              G()[c*dim() + delta] = detJ*sum; 
            }
        }
    }
  addFlops(flops);
}



void QuadratureIntegral
::createOneFormTransformationMatrix(const CellJacobianBatch& J,  
                                    const Array<int>& alpha) const 
{
  TimeMonitor timer(quadTransCreationTimer());

  TEST_FOR_EXCEPTION(J.cellDim() != dim(), InternalError,
                     "Inconsistency between Jacobian dimension " << J.cellDim()
                     << " and cell dimension " << dim() 
                     << " in QuadratureIntegral::createOneFormTransformationMatrix()");
  G().resize(J.numCells() * J.cellDim());

  int flops = J.numCells() * dim() * (1 + alpha.size());

  for (int c=0; c<J.numCells(); c++)
    {

      Array<double> invJ;
      J.getInvJ(c, invJ);
      double detJ = fabs(J.detJ()[c]);
      for (int gamma=0; gamma<dim(); gamma++)
        {
          double sum = 0.0;
          for (int t=0; t<alpha.size(); t++)
            {
              sum += invJ[alpha[t] + gamma*dim()];
            }
          G()[c*dim() + gamma] = detJ*sum; 
        }
    }
  
  addFlops(flops);
}
