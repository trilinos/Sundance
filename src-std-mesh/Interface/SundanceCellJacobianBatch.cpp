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
  void dgetrs_(const char* trans, const int* N, const int* NRHS, 
               const double* A, const int* lda, 
               const int* iPiv, double* B, const int* ldb, int* info);

  /* LAPACK factorization */
  void dgetrf_(const int* M, const int* N, double* A, const int* lda, 
               const int* iPiv, int* info);
};


CellJacobianBatch::CellJacobianBatch()
  : spatialDim_(0), cellDim_(0), 
    jSize_(0), numCells_(0), numQuad_(0), iPiv_(), J_(), detJ_(), invJ_(),
    isFactored_(false), hasInverses_(false)
{;}

void CellJacobianBatch::resize(int numCells, int numQuad, 
                               int spatialDim, int cellDim)
{
  spatialDim_ = spatialDim;
  cellDim_ = cellDim;
  if (spatialDim_ == cellDim_)
    {
      jSize_ = spatialDim_*spatialDim_;
     }
  
  numCells_ = numCells;
  numQuad_ = numQuad;
  iPiv_.resize(spatialDim_*numCells_*numQuad_);
  J_.resize(spatialDim_*spatialDim_*numCells_*numQuad_);
  detJ_.resize(numCells_*numQuad_);
  isFactored_ = false;
  hasInverses_ = false;
}

void CellJacobianBatch::resize(int numCells, int spatialDim, int cellDim)
{
  spatialDim_ = spatialDim;
  cellDim_ = cellDim;
  if (spatialDim_ == cellDim_)
    {
      jSize_ = spatialDim_*spatialDim_;
    }

  numCells_ = numCells;
  numQuad_ = 1;
  iPiv_.resize(spatialDim_*numCells_);
  J_.resize(spatialDim_*spatialDim_*numCells_);
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

  TEST_FOR_EXCEPTION(spatialDim_ != cellDim_, InternalError,
                     "Attempting to factor the Jacobian of a cell "
                     "that is not of maximal dimension");
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
          int* iPiv = &(iPiv_[(q + cell*numQuad_)*spatialDim_]);
  
          /* fortran junk */
          int lda = spatialDim_; // leading dimension of J
          
          int info = 0; // error return flag, will be zero if successful. 
          
          /* Factor J */
          ::dgetrf_( &spatialDim_,  &spatialDim_, jFactPtr, &lda, iPiv, &info);
          
          TEST_FOR_EXCEPTION(info != 0, RuntimeError,
                             "CellJacobianBatch::setJacobian(): factoring failed");

          /* the determinant is the product of the diagonal elements 
           * the upper triangular factor of the factored Jacobian */
          double detJ = 1.0;
          for (int i=0; i<spatialDim_; i++)
            {
              detJ *= jFactPtr[i + spatialDim_*i];
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

  invJ_.resize(spatialDim_*spatialDim_*numQuad_*numCells_);

  if (!isFactored_) factor();
  
  for (int cell=0; cell<numCells_; cell++)
    {
      for (int q=0; q<numQuad_; q++)
        {
          int start = (cell*numQuad_ + q)*jSize_;

          /* pointer to start of J for this cell */
          double* jFactPtr = &(J_[start]);
          double* invJPtr = &(invJ_[start]);
          int* iPiv = &(iPiv_[(q + cell*numQuad_)*spatialDim_]);
  
          /* fortran junk */
          int lda = spatialDim_; // leading dimension of J
          
          int info = 0; // error return flag, will be zero if successful. 
          
          /* fill the inverse of J with the identity */
          for (int i=0; i<spatialDim_; i++)
            {
              for (int j=0; j<spatialDim_; j++)
                {
                  if (i==j) invJPtr[i*spatialDim_+j] = 1.0;
                  else invJPtr[i*spatialDim_+j] = 0.0;
                }
            }

          ::dgetrs_("N",  &spatialDim_,  &spatialDim_, jFactPtr, 
                     &spatialDim_, iPiv, invJPtr,  &spatialDim_, &info);
          
          TEST_FOR_EXCEPTION(info != 0, RuntimeError,
                             "CellJacobianBatch::setJacobian(): inversion failed");
        }
    }
}

void CellJacobianBatch::applyInvJ(int cell, int q, 
                                  double* rhs, int nRhs, bool trans) const 
{
  if (!isFactored_) factor();

  double* jFactPtr = &(J_[(cell*numQuad_ + q)*spatialDim_*spatialDim_]);
  int* iPiv = &(iPiv_[(q + cell*numQuad_)*spatialDim_]);

  int info = 0; // error return flag, will be zero if successful. 
  
  if (trans)
    {
      ::dgetrs_("T",  &spatialDim_, &nRhs, jFactPtr,  &spatialDim_, 
                iPiv, rhs,  &spatialDim_, &info);
    }
  else
    {
      ::dgetrs_("N",  &spatialDim_, &nRhs, jFactPtr,  &spatialDim_, 
                iPiv, rhs,  &spatialDim_, &info);
    }
          
  TEST_FOR_EXCEPTION(info != 0, RuntimeError,
                     "CellJacobianBatch::applyInvJ(): backsolve failed");
}

void CellJacobianBatch::getInvJ(int cell, int quad, Array<double>& invJ) const 
{
  if (!hasInverses_) computeInverses();
  
  int start = (cell*numQuad_ + quad)*jSize_;
  
  invJ.resize(spatialDim_*spatialDim_);

  for (int col=0; col<spatialDim_; col++)
    {
      for (int row=0; row<spatialDim_; row++) 
        {
          invJ[col + spatialDim_*row] = invJ_[start + col + spatialDim_*row];
        }
    }
}



