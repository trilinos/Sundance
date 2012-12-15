/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_IFPACKICCOPERATOR_HPP
#define PLAYA_IFPACKICCOPERATOR_HPP


#include "PlayaEpetraMatrix.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "Ifpack_ICT.h"


namespace Playa
{
/**
 * This is the operator representation of the inverse upper triangular factor 
 * \f${\tilde R}^{-1}\f$  
 * in the incomplete Cholesky factorization of a matrix \f$ A \approx {\tilde R}{\tilde R}^T\f$. 
 *
 * @author Kimberly Kennedy, with some small modifications by Kevin Long
 */
class IfpackICCOperator : 
  public LinearOpWithSpaces<double>,
  public Printable
{
public:
  /** Construct the operator. During construction, the matrix A will be
   * approximately factored using the options given in the arguments. 
   *
   * \param A The matrix to be factored
   * \param fillLevels Specifies the number of on-processor 
   * nonzeros allowed in each
   * row of the factorization, given as a ratio of the 
   * allowed nnz per row in the factor to nnz per row in A. 
   * \param overlapFill Fill allowed for off-processor elements
   * \param dropTol New elements are dropped if |R_ij| < dropTol*R_jj.
   * \param relaxationValue Fraction of dropped element mass to be added to 
   * diagonal.
   * \param relativeThreshold Fraction of diagonal element to be added to diagonal
   * \param absoluteThreshold Amount to be added to each diagonal element
   */
  IfpackICCOperator(const EpetraMatrix* A,
		    int fillLevels,
		    int overlapFill,
		    double dropTol,
		    double relaxationValue,
		    double relativeThreshold,
		    double absoluteThreshold);

  /** 
   * Apply the operator. 
   */
  virtual void apply(
		     Teuchos::ETransp applyType,
		     const Vector<double>& in,
		     Vector<double> out) const ;
  
  /** \name Diagnostic output */
  //@{
  /** Print the matrix */
  virtual void print(std::ostream& os) const ;
  //@}

  /** */
  std::ostream& describe(
			 std::ostream                         &out
			 ,const Teuchos::EVerbosityLevel      verbLevel
			 ,const std::string                   leadingIndent
			 , const std::string                   indentSpacer
			 ) const 
  {
    out << leadingIndent << indentSpacer << this->description() << std::endl;
    return out;
  }
  /** */
  std::string description() const ;
  //@}



private:
  RCP<Ifpack_ICT> precond_;

};


}

#endif 
