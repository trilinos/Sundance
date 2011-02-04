/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_INCREMENTALLYCONFIGURABLEMATRIXFACTORY_HPP
#define PLAYA_INCREMENTALLYCONFIGURABLEMATRIXFACTORY_HPP

#include "PlayaDefs.hpp"

namespace Playa
{
  /** 
   * Class IncrementallyConfigurableMatrixFactory provides an abstract 
   * interface for row-at-a-time configuration of matrix factories.
   */
  class IncrementallyConfigurableMatrixFactory
  {
  public:
    /** Virtual dtor */
    virtual ~IncrementallyConfigurableMatrixFactory(){;}

    /** Initialize a set of nonzero elements in the matrix's graph.
     * @param globalRowIndex the global index of the row to which these
     * elements belong.
     * @param nElemsToInsert the number of elements being inserted in this
     * step
     * @param globalColumnIndices array of column indices. Must 
     * be nElemsToInsert in length. 
     */
    virtual void initializeNonzerosInRow(int globalRowIndex,
                                         int nElemsToInsert,
                                         const int* globalColumnIndices) = 0 ;

    /** 
     * Initialize nonzeros in a batch of rows. 
     */
    virtual void initializeNonzeroBatch(int numRows, 
                                        int rowBlockSize,
                                        const int* globalRowIndices,
                                        int numColumnsPerRow,
                                        const int* globalColumnIndices,
                                        const int* skipRow);

    /** Finalize values of the matrix. This is a hook for any
     * implementation-dependent steps that must be done after
     * loading of elements. */
    virtual void finalize() = 0 ;

  private:
    
    
  };

  /* Default implementation of initializeElementBatch */
  inline void IncrementallyConfigurableMatrixFactory
  ::initializeNonzeroBatch(int numRows, 
                           int rowBlockSize,
                           const int* globalRowIndices,
                           int numColumnsPerRow,
                           const int* globalColumnIndices,
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
            initializeNonzerosInRow(globalRowIndices[row], 
                                    numColumnsPerRow, cols);
          }
      }
  }
}

#endif
