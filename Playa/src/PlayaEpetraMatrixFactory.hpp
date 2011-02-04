/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_EPETRAMATRIXFACTORY_HPP
#define PLAYA_EPETRAMATRIXFACTORY_HPP

#include "PlayaEpetraVectorSpace.hpp"
#include "PlayaIncrementallyConfigurableMatrixFactory.hpp"
#include "PlayaCollectivelyConfigurableMatrixFactory.hpp"
#include "PlayaMatrixFactory.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "Epetra_CrsGraph.h"

namespace Playa
{
using namespace Teuchos;
  

/** */
class EpetraMatrixFactory : public MatrixFactory<double>,
                            public IncrementallyConfigurableMatrixFactory,
                            public CollectivelyConfigurableMatrixFactory
{
public:

  /** Construct an uninitialized EpetraMatrixFactory */
  EpetraMatrixFactory(const RCP<const EpetraVectorSpace>& domain,
    const RCP<const EpetraVectorSpace>& range);

  /** */
  const RCP<const EpetraVectorSpace>& epRange() const {return range_;}

  /** */
  const RCP<const EpetraVectorSpace>& epDomain() const {return domain_;}


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
    const int* globalColumnIndices) ;

  /** 
   * Initialize nonzeros in a batch of rows. 
   */
  virtual void initializeNonzeroBatch(int numRows, 
    int rowBlockSize,
    const int* globalRowIndices,
    int numColumnsPerRow,
    const int* globalColumnIndices,
    const int* skipRow);

  /** Configure all rows at once */
  virtual void configure(int lowestRow,
    const std::vector<int>& rowPtrs,
    const std::vector<int>& nnzPerRow,
    const std::vector<int>& data);

  /** */
  void finalize();

  /** */
  const Epetra_CrsGraph& graph() const ;

  /** */
  virtual LinearOperator<double> createMatrix() const ;

protected:

private:

  /** */
  RCP<Epetra_CrsGraph> graph_;

  /** */
  RCP<const EpetraVectorSpace> range_;

  /** */
  RCP<const EpetraVectorSpace> domain_;
};
}

#endif
