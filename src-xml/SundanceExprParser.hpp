#ifndef SUNDANCE_EXPRPARSER_H
#define SUNDANCE_EXPRPARSER_H

#include "SundanceDefs.hpp"
#include "SundanceToken.hpp"
#include "SundanceSet.hpp"
#include "SundanceExprScanner.hpp"
#include "Teuchos_XMLObject.hpp"


namespace SundanceXML
{
  class ExprParser
  {
  public:

    static Teuchos::XMLObject evaluate(const string& line);
    static Teuchos::XMLObject evaluateAssignment(const string& line);


    static void assignLevel(const ExprScanner& scanner,
                            Teuchos::XMLObject& result);
    static void addLevel(const ExprScanner& scanner,
                         Teuchos::XMLObject& result);
    static void multLevel(const ExprScanner& scanner,
                          Teuchos::XMLObject& result);
    static void unaryLevel(const ExprScanner& scanner,
                           Teuchos::XMLObject& result);
    static void parenLevel(const ExprScanner& scanner,
                           Teuchos::XMLObject& result);
    static void primitiveLevel(const ExprScanner& scanner,
                               Teuchos::XMLObject& result);

    static void parseFunction(const string& funcName,
                              const Teuchos::XMLObject& arg,
                              Teuchos::XMLObject& result);

    static void parseExpr(const string& funcName,
                          const Teuchos::XMLObject& arg,
                          Teuchos::XMLObject& result); 

    static void parseMesh(const string& funcName,
                          const Teuchos::XMLObject& arg,
                          Teuchos::XMLObject& result); 

    static void parseBasis(const string& funcName,
                           const Teuchos::XMLObject& arg,
                           Teuchos::XMLObject& result);

    static void parseQuad(const string& funcName,
                          const Teuchos::XMLObject& arg,
                          Teuchos::XMLObject& result);

    static void parseDiscreteSpace(const string& funcName,
                                   const Teuchos::XMLObject& arg,
                                   Teuchos::XMLObject& result); 

    static void parseLinearProblem(const string& funcName,
                                   const Teuchos::XMLObject& arg,
                                   Teuchos::XMLObject& result);

    static void parseNonlinearProblem(const string& funcName,
                                      const Teuchos::XMLObject& arg,
                                      Teuchos::XMLObject& result);

    static SundanceUtils::Set<string>& exprTypes() ;

    static SundanceUtils::Set<string>& basisTypes() ;

    static SundanceUtils::Set<string>& quadTypes() ;

    static SundanceUtils::Set<string>& meshTypes() ;
  private:

  };

}
#endif
