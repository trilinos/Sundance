#ifndef SUNDANCE_INTERPRETER_H
#define SUNDANCE_INTERPRETER_H

#include "Sundance.hpp"

namespace SundanceXML
{
  class Interpreter
  {
  public:

    void interpret(const XMLObject& xml);

    const Mesh& mesh(const string& name) const ;

    const Expr& expr(const string& name) const ;

    const QuadratureFamily& quad(const string& name) const ;

    const BasisFamily& basis(const string& name) const ;

    const DiscreteSpace& expr(const string& name) const ;

    const LinearProblem& linearProblem(const string& name) const ;
  private:
    
    SundanceUtils::Map<string, Mesh> meshes_;

    SundanceUtils::Map<string, Expr> exprs_;

    SundanceUtils::Map<
  };

}
#endif
