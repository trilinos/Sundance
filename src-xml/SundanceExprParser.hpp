#ifndef SUNDANCE_EXPRPARSER_H
#define SUNDANCE_EXPRPARSER_H

#include "SundanceDefs.hpp"
#include "SundanceToken.hpp"
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



  private:

  };

}
#endif
