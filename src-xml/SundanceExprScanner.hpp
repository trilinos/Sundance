#ifndef SUNDANCE_EXPRSCANNER_H
#define SUNDANCE_EXPRSCANNER_H

#include "SundanceDefs.hpp"
#include "SundanceToken.hpp"
#include "Teuchos_Array.hpp"

namespace SundanceXML
{
  class ExprScanner
    {
    public:
      ExprScanner(const string& line);

      const Token& pop() const ;
      const Token& peek() const ;
      const Token& peekAtNext() const ;
      bool hasMore() const ;

      string showError() const ;

    private:
      bool isLogical(char c) const ;
      bool isWhite(char c) const ;
      bool isDelimiter(char c) const ;
      bool isQuote(char c) const ;
      Token getToken(const string& line, int& place) const ;

      Teuchos::Array<Token> tok_;
      mutable int ptr_;
    };


}
#endif
