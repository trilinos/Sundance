/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_ROWACCESSIBLEOP_HPP
#define PLAYA_ROWACCESSIBLEOP_HPP

#include "PlayaDefs.hpp"
#include "Teuchos_Array.hpp"

namespace Playa
{
  /** 
   * Base interface for operators for which a row may be extracted.
   */
  template <class Scalar>
  class RowAccessibleOp 
  {
  public:
    /** Virtual dtor */
    virtual ~RowAccessibleOp(){;}

    /** 
     * Get the non-zero values in the row-th row.
     * @param row the index of the row
     * @param indices the column indices of the non-zero values in row row
     * @param values the non-zero values corresponding to the indices in indices
     */
    virtual void getRow(const int& row, 
			Teuchos::Array<int>& indices, 
			Teuchos::Array<Scalar>& values) const = 0;

  private:
    
    
  };
}

#endif
