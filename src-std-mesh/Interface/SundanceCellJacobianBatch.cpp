/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellJacobianBatch.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceStdMesh::Internal;
using namespace SundanceStdMesh;
using namespace Teuchos;
using namespace SundanceUtils;
using namespace TSFExtended;

/* declare LAPACK subroutines */
extern "C"
{
  /* LAPACK backsolve on a factored system */
  void dgetrs_(char* trans, int* N, int* NRHS, double* A, int* lda, 
               int* iPiv, double* B, int* ldb, int* info);

  /* LAPACK factorization */
  void dgetrf_(int* M, int* N, double* A, int* lda, int* iPiv, int* info);
};


CellJacobianBatch::CellJacobianBatch()
  : dim_(0), jSize_(0), numCells_(0), numQuad_(0), iPiv_(), J_(), detJ_(), invJ_(),
    isFactored_(false), hasInverses_(false)
{;}

void CellJacobianBatch::resize(int numCells, int numQuad, int dim)
{
  dim_ = dim;
  jSize_ = dim_*dim_;
  numCells_ = numCells;
  numQuad_ = numQuad;
  iPiv_.resize(dim_*numCells_*numQuad_);
  J_.resize(dim_*dim_*numCells_*numQuad_);
  detJ_.resize(numCells_*numQuad_);
  isFactored_ = false;
  hasInverses_ = false;
}

void CellJacobianBatch::resize(int numCells, int dim)
{
  dim_ = dim;
  jSize_ = dim_*dim_;
  numCells_ = numCells;
  numQuad_ = 1;
  iPiv_.resize(dim_*numCells_);
  J_.resize(dim_*dim_*numCells_);
  detJ_.resize(numCells_);
  isFactored_ = false;
  hasInverses_ = false;
}

void CellJacobianBatch::factor() const 
{
  if (isFactored_) return;
  /* We're given the Jacobian, and we want to factor it and compute its determinant. 
   * We factor it using the LAPACK routine dgetrf(), after which J is replaced
   * by its LU factorization. The determinant of J is obtained by taking the
   * project of the diagonal elements of U. 
   */
 
  Tabs tabs;
  SUNDANCE_OUT(verbosity() > VerbMedium,
               tabs << "factoring Jacobians");
  
  for (int cell=0; cell<numCells_; cell++)
    {
      for (int q=0; q<numQuad_; q++)
        {
          int start = (cell*numQuad_ + q)*jSize_;

          /* pointer to start of J for this cell */
          double* jFactPtr = &(J_[start]);
          int* iPiv = &(iPiv_[(q + cell*numQuad_)*dim_]);
  
          /* fortran junk */
          int lda = dim_; // leading dimension of J
          
          int info = 0; // error return flag, will be zero if successful. 
          
          /* Factor J */
          ::dgetrf_((int*) &dim_, (int*) &dim_, jFactPtr, &lda, iPiv, &info);
          
          TEST_FOR_EXCEPTION(info != 0, RuntimeError,
                             "CellJacobianBatch::setJacobian(): factoring failed");

          /* the determinant is the product of the diagonal elements 
           * the upper triangular factor of the factored Jacobian */
          double detJ = 1.0;
          for (int i=0; i<dim_; i++)
            {
              detJ *= jFactPtr[i + dim_*i];
            }
          detJ_[cell*numQuad_ + q] = detJ;
        }
    }

  isFactored_ = true;
}

void CellJacobianBatch::computeInverses() const 
{
  if (hasInverses_) return;

  Tabs tabs;
  SUNDANCE_OUT(verbosity() > VerbMedium,
               tabs << "inverting Jacobians");

  invJ_.resize(dim_*dim_*numQuad_*numCells_);

  if (!isFactored_) factor();
  
  for (int cell=0; cell<numCells_; cell++)
    {
      for (int q=0; q<numQuad_; q++)
        {
          int start = (cell*numQuad_ + q)*jSize_;

          /* pointer to start of J for this cell */
          double* jFactPtr = &(J_[start]);
          double* invJPtr = &(invJ_[start]);
          int* iPiv = &(iPiv_[(q + cell*numQuad_)*dim_]);
  
          /* fortran junk */
          int lda = dim_; // leading dimension of J
          
          int info = 0; // error return flag, will be zero if successful. 
          
          /* fill the inverse of J with the identity */
          for (int i=0; i<dim_; i++)
            {
              for (int j=0; j<dim_; j++)
                {
                  if (i==j) invJPtr[i*dim_+j] = 1.0;
                  else invJPtr[i*dim_+j] = 0.0;
                }
            }

          ::dgetrs_("N", (int*) &dim_, (int*) &dim_, jFactPtr, 
                    (int*) &dim_, iPiv, invJPtr, (int*) &dim_, &info);
          
          TEST_FOR_EXCEPTION(info != 0, RuntimeError,
                             "CellJacobianBatch::setJacobian(): inversion failed");
        }
    }
}

void CellJacobianBatch::applyInvJ(int cell, int q, 
                                  double* rhs, int nRhs, bool trans) const 
{
  if (!isFactored_) factor();

  double* jFactPtr = &(J_[(cell*numQuad_ + q)*dim_*dim_]);
  int* iPiv = &(iPiv_[(q + cell*numQuad_)*dim_]);

  int info = 0; // error return flag, will be zero if successful. 
  
  if (trans)
    {
      ::dgetrs_("T", (int*) &dim_, &nRhs, jFactPtr, (int*) &dim_, iPiv, rhs, (int*) &dim_, &info);
    }
  else
    {
      ::dgetrs_("N", (int*) &dim_, &nRhs, jFactPtr, (int*) &dim_, iPiv, rhs, (int*) &dim_, &info);
    }
          
  TEST_FOR_EXCEPTION(info != 0, RuntimeError,
                     "CellJacobianBatch::applyInvJ(): backsolve failed");
}

void CellJacobianBatch::getInvJ(int cell, int quad, Array<double>& invJ) const 
{
  if (!hasInverses_) computeInverses();
  
  int start = (cell*numQuad_ + quad)*jSize_;
  
  invJ.resize(dim_*dim_);

  for (int col=0; col<dim_; col++)
    {
      for (int row=0; row<dim_; row++) 
        {
          invJ[col + dim_*row] = invJ_[start + col + dim_*row];
        }
    }
}




