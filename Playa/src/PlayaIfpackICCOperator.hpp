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
 *
 */
class IfpackICCOperator : 
  public LinearOpWithSpaces<double>,
  public Printable
{
public:
  /** */
  IfpackICCOperator(const EpetraMatrix* A,
    int fillLevels,
    int overlapFill,
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
