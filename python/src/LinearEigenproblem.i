// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceLinearEigenproblem.hpp"
#include "PySundanceLinearSolver.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "std_complex.i"
%include "exception.i"

namespace Sundance
{

class Eigensolution
{
public:
  /** */
  bool isValid() const ;

  /** */
  int numEigenfunctions() const ;

  /** */
  const Sundance::Expr& eigenfunction(int i) const ;

  /** */
  const std::complex<double>& eigenvalue(int i) const 
    {return eigenvalues_[i];}

};
}

%rename(ComplexArray) CompArray;
%rename(VectorArray) VecArray;
%template(VecArray) Teuchos::Array<TSFExtended::Vector<double> >;
%template(CompArray) Teuchos::Array<std::complex<double> >;

namespace TSFExtended
{
%rename(Eigensolver) EigenSol;

  template <class Scalar>
  class Eigensolver
  {
  public:
    Eigensolver();

    /** */
    void solve(
      const LinearOperator<Scalar>& K,
      const LinearOperator<Scalar>& M,
      Teuchos::Array<Vector<Scalar> >& ev,
      Teuchos::Array<std::complex<Scalar> >& ew) const ;

    /** */
    void solve(
    const LinearOperator<Scalar>& K,
    Teuchos::Array<Vector<Scalar> >& ev,
    Teuchos::Array<std::complex<Scalar> >& ew) const ;

  };

%template(EigenSol) Eigensolver<double>;
}



 
namespace Sundance
{ 
class LinearEigenproblem
{
public:
  /** */
  LinearEigenproblem(const Sundance::Mesh& mesh, 
    const Sundance::Expr& eqn,
    const Sundance::Expr& v, 
    const Sundance::Expr& u,
    const TSFExtended::VectorType<double>& vecType);
  /** */
  LinearEigenproblem(
    const Sundance::Mesh& mesh, const Sundance::Expr& eqn,
    const Sundance::Expr& v, const Sundance::Expr& u,
    const TSFExtended::VectorType<double>& vecType,
    bool lumpMass) ;
  /** */
  LinearEigenproblem(
    const Sundance::Mesh& mesh, const Sundance::Expr& eqn,
    const Sundance::Expr& massExpr,
    const Sundance::Expr& v, const Sundance::Expr& u,
    const TSFExtended::VectorType<double>& vecType,
    bool lumpMass) ;
    
  /** */
  Eigensolution solve(const TSFExtended::Eigensolver<double>& solver) const ;

  /** */
  TSFExtended::LinearOperator<double> getK() const ;
      

  /** */
  TSFExtended::LinearOperator<double> getM() const ;
    

};


}
