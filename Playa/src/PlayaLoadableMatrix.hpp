/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_LOADABLEMATRIX_HPP
#define PLAYA_LOADABLEMATRIX_HPP

#include "PlayaDefs.hpp"

namespace Playa
{
  /** 
   * Class LoadableMatrix provides an abstract interface for 
   * loading elements into a matrix.
   */
  template <class Scalar>
  class LoadableMatrix 
  {
  public:
    /** Virtual dtor */
    virtual ~LoadableMatrix(){;}

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
                          const Scalar* elementValues) = 0 ;

    /** Set all elements to zero, preserving the existing structure */
    virtual void zero() = 0 ;

    /** 
     * Add to a batch of elements
     */
    virtual void addToElementBatch(int numRows, 
                                   int rowBlockSize,
                                   const int* globalRowIndices,
                                   int numColumnsPerRow,
                                   const int* globalColumnIndices,
                                   const Scalar* values,
                                   const int* skipRow);


  };


  /* Default implementation of addElementBatch */
  template <class Scalar>
  void LoadableMatrix<Scalar>::addToElementBatch(int numRows, 
                                                 int rowBlockSize,
                                                 const int* globalRowIndices,
                                                 int numColumnsPerRow,
                                                 const int* globalColumnIndices,
                                                 const Scalar* values,
                                                 const int* skipRow)
  {
    int numRowBlocks = numRows/rowBlockSize;
    int row = 0;

    for (int rb=0; rb<numRowBlocks; rb++)
      {
        const int* cols = globalColumnIndices + rb*numColumnsPerRow;
        for (int r=0; r<rowBlockSize; r++, row++)
          {
            if (skipRow[row]) continue;
            const double* rowVals = values + row*numColumnsPerRow;
            addToRow(globalRowIndices[row], numColumnsPerRow,
                     cols, rowVals);
          }
      }
  }
}

#endif
