/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_ILUFACTORIZABLEOP_HPP
#define PLAYA_ILUFACTORIZABLEOP_HPP

#include "PlayaDefs.hpp"


namespace Playa
{
template <class Scalar> class Preconditioner;

/** */
enum LeftOrRight {Left, Right};

/** 
 * Base interface for operators for which incomplete LU factorizations
 * can be obtained. 
 */
template <class Scalar>
class ILUFactorizableOp
{
public:
  /** Virtual dtor */
  virtual ~ILUFactorizableOp(){;}


  /** \name incomplete factorization preconditioning interface */
  //@{
  /** create an incomplete factorization. 
   * @param fillLevels number of levels of fill on the local processor
   * @param overlapFill number of levels of fill on remote processors
   * @param relaxationValue fraction of dropped values to be added to the
   * diagonal
   * @param relativeThreshold relative diagonal perutrbation
   * @param absoluteThreshold absolute diagonal perturbation
   * @param leftOrRight whether this preconditioner is to be applied
   * from the left or right 
   * @param rtn newly created preconditioner, returned 
   * by reference argument.
   */
  virtual void getILUKPreconditioner(int fillLevels,
    int overlapFill,
    double relaxationValue,
    double relativeThreshold,
    double absoluteThreshold,
    LeftOrRight leftOrRight,
    Preconditioner<Scalar>& rtn) const=0;
  //@}
     
      
private:
};
}


#endif
