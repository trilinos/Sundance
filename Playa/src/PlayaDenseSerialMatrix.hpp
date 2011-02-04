/* @HEADER@ */
//
 /* @HEADER@ */



#ifndef PLAYA_DENSE_SERIAL_MATRIX_H
#define PLAYA_DENSE_SERIAL_MATRIX_H

#include "PlayaDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "PlayaPrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "PlayaSerialVectorSpace.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaSolverState.hpp"

namespace Playa
{
using namespace Teuchos;

template <class T> class LinearOperator;

/**
 * Linear operator implemented as a dense matrix.
 */

class DenseSerialMatrix : public LinearOpWithSpaces<double>,
                          public LoadableMatrix<double>,
                          public Playa::Printable
{
public:
  /** Construct with domain and range spaces, which should be
   * DenseSerialVectorSpace objects */
  DenseSerialMatrix(
    const VectorSpace<double>& domain,
    const VectorSpace<double>& range);

  /** Virtual dtor */
  ~DenseSerialMatrix(){;}

  /** 
   * Apply either the operator or its transpose
   */
  virtual void apply(Teuchos::ETransp transApplyType,
    const Vector<double>& in,
    Vector<double> out) const ;

  /** Insert a set of elements in a row, adding to any previously
   * existing values.  The nonzero structure of the matrix must have
   * been determined at construction time. 
   *
   * @param globalRowIndex the global index of the row to which these
   * elements belong.
   * @param nElemsToInsert the number of elements being inserted in this
   * step
   * @param globalColumnIndices array of column indices. Must 
   * be nElemsToInsert in length. 
   * @param elements array of element values. Must be nElemsToInsert in
   * length
   */
  virtual void addToRow(int globalRowIndex,
    int nElemsToInsert,
    const int* globalColumnIndices,
    const double* elementValues) ;

  /** Set all elements to zero, preserving the existing structure */
  virtual void zero() ;



  /** write to a stream */
  void print(std::ostream& os) const ;

  /** */
  const double * const dataPtr() const {return &(data_[0]);}

  /** */
  double* dataPtr() {return &(data_[0]);}

  /** */
  int numRows() const {return nRows_;}

  /** */
  int numCols() const {return nCols_;}

  /** */
  void setRow(int row, const Array<double>& rowVals);


private:

  int nRows_;
  int nCols_;
  Array<double> data_;
};


/** \relates DenseSerialMatrix */
void denseSVD(const LinearOperator<double>& A,
  LinearOperator<double>& U,  
  Vector<double>& Sigma,
  LinearOperator<double>& Vt);

/** \relates DenseSerialMatrix */
SolverState<double> denseSolve(const LinearOperator<double>& A,
  const Vector<double>& b,
  Vector<double>& x);

}

#endif
