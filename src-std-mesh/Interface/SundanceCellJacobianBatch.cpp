/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellJacobianBatch.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdMesh::Internal;
using namespace SundanceStdMesh;
using namespace Teuchos;
using namespace SundanceUtils;

/* declare dgesv(), the LAPACK linear system solver */
extern "C"
{
  void dgesv_(int *n, int *nrhs, double *a, int *lda, 
              int *ipiv, double *b, int *ldb, int *info);
};


CellJacobianBatch::CellJacobianBatch()
  : dim_(0), numCells_(0), numQuad_(0), invJ_(), detJ_()
{;}

void CellJacobianBatch::resize(int numCells, int numQuad, int dim)
{
  dim_ = dim;
  numCells_ = numCells;
  numQuad_ = numQuad;
  invJ_.resize(dim_*dim_*numCells_*numQuad_);
  detJ_.resize(numCells_*numQuad_);
}

void CellJacobianBatch::resize(int numCells, int dim)
{
  dim_ = dim;
  numCells_ = numCells;
  numQuad_ = 1;
  invJ_.resize(dim_*dim_*numCells_);
  detJ_.resize(numCells_);
}

void CellJacobianBatch::setJacobian(int cell, int q, Array<double>& jVals)
{
  /* We're given the Jacobian, and we want to compute its inverse and determinant 
   * We can get the inverse by solving the linear system J J^-1 = I using the
   * LAPACK routine dgesv. In calling dgesv, the input Jacobian is overwritten
   * by its LU factorization. The determinant of J is obtained by taking the
   * project of the diagonal elements of U */
 
  if (verbose())
    {
      cerr << "setting Jacobian for cell=" << cell 
           << " quad pt = " << q << endl;
      viewMatrix(cerr, "J", &(jVals[0]));
    }

  /* keep around a pivot vector */
  static Array<int> iPiv(3);

  int jSize = dim_*dim_;
  int start = (cell*numQuad_ + q)*jSize;

  /* pointer to start of input jacobian. Upon calling dgesv, this 
   * will be overwritten with the LU factorization of J */
  double* jPtr  = &(jVals[0]);
  /* pointer to start of invJ */
  double* jInvPtr = &(invJ_[start]);
  double* jp = jInvPtr;
  
  /* initialize J^-1 as the identity, to use as the RHS in the call to dgesv */
  for (int i=0; i<dim_; i++)
    {
      for (int j=0; j<dim_; j++, jp++)
        {
          if (i==j) *jp = 1.0;
          else *jp = 0.0;
        }
    }

  /* fortran junk */
  int lda = dim_; // leading dimension of J
  int ldb = dim_; // leading dimension of RHS 

  int info = 0; // error return flag, will be zero if successful. 

  /* solve J J^-1 = I for J^-1 */
  ::dgesv_(&dim_, &dim_, jPtr, &lda, &(iPiv[0]), jInvPtr, &ldb, &info);

  TEST_FOR_EXCEPTION(info != 0, RuntimeError,
                     "CellJacobianBatch::setJacobian(): dgesv failed");

  /* the determinant is the product of the diagonal elements 
   * the upper triangular factor of the factored Jacobian */
  double detJ = 1.0;
  for (int i=0; i<dim_; i++)
    {
      detJ *= jPtr[i + dim_*i];
    }
  detJ_[cell*numQuad_ + q] = detJ;

  if (verbose())
    {
      viewMatrix(cerr, "J^-1", jInvPtr);
      cerr << "det J = " << detJ << endl;
    }
}


void CellJacobianBatch::viewMatrix(ostream& os, 
                               const string& header, const double* vals) const
{
  os << header << " = {";
  int k = 0;
  for (int i=0; i<dim_; i++)
    {
      os << "{";
      for (int j=0; j<dim_; j++, k++)
        {
          os << vals[k];
          if (j < dim_-1) os << ", ";
        }
      os << "}";
    }
  os << "}" << endl;
}



