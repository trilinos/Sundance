#ifndef SUNDANCE_INTERPRETER_H
#define SUNDANCE_INTERPRETER_H

#include "Sundance.hpp"

#define LOOKUP(typeName, mapName, objName) \
typeName rtn; \
TEST_FOR_EXCEPTION(!mapName.containsKey(objName), RuntimeError, \
"object [" << objName << " not found in map " << mapName.toString()); \
return mapName.get(objName);

namespace SundanceXML
{
  class Interpreter
  {
  public:

    void interpret(const XMLObject& xml);

    const Mesh& mesh(const string& name) const {LOOKUP(Mesh, mesh_, name);}

    const Expr& expr(const string& name) const {LOOKUP(Expr, expr_, name);}

    const Expr& unknown(const string& name) const {LOOKUP(Expr, unk_, name);}

    const Expr& test(const string& name) const {LOOKUP(Expr, test_, name);}

    const QuadratureFamily& quad(const string& name) const ;

    const BasisFamily& basis(const string& name) const ;

    const DiscreteSpace& expr(const string& name) const ;

    const LinearProblem& linearProblem(const string& name) const ;
  private:
    
    SundanceUtils::Map<string, Mesh> mesh_;

    SundanceUtils::Map<string, Expr> expr_;

    SundanceUtils::Map<string, Expr> unk_;

    SundanceUtils::Map<string, Expr> test_;

    SundanceUtils::Map<string, BasisFamily> basis_;

    SundanceUtils::Map<string, QuadratureFamily> basis_;

    SundanceUtils::Map<string, CellFilter> filter_;

  };

}
#endif
