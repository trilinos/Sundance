/* @HEADER@ */
//   
 /* @HEADER@ */


#ifndef MATRIXLAPLACIAN1D_HPP
#define MATRIXLAPLACIAN1D_HPP

#include "PlayaOperatorBuilder.hpp"

using namespace Playa;
using namespace Teuchos;


namespace Playa
{

/** */
class MatrixLaplacian1D : public OperatorBuilder<double>
{
public:
    
  /** */
  MatrixLaplacian1D(int nLocal, const VectorType<double>& vecType);

  /** */
  LinearOperator<double> getOp() const {return op_;}

private:
  /** */
  void init(int nLocal, const VectorType<double>& vecType);
  
  LinearOperator<double> op_;
};

}

#endif
