#ifndef SUNDANCE_TOKEN_H
#define SUNDANCE_TOKEN_H

#include "SundanceDefs.hpp"

namespace SundanceXML
{
  enum TokenType {TT_Unknown, TT_Assign, TT_OpenParen, TT_CloseParen,
                  TT_OpenBracket, TT_CloseBracket, TT_OpenBrace,
                  TT_CloseBrace, TT_Comma, TT_Plus, TT_Minus, TT_Times,
                  TT_Divide, TT_Constant, TT_Name, TT_Equal,
                  TT_LessThan, TT_GreaterThan, TT_LessThanOrEqual,
                  TT_GreaterThanOrEqual, TT_And, TT_Or,
                  TT_Not, TT_NotEqual, TT_QuotedString, TT_Dot, TT_End};

  class Token
    {
    public:
      Token() : tok_(), type_(TT_Unknown) {;}
      Token(const string& tok);

      bool isUnknown() const {return type_==TT_Unknown;}

      bool isAssign() const {return type_==TT_Assign;}
      bool isOpenParen() const {return type_==TT_OpenParen;}
      bool isCloseParen() const {return type_==TT_CloseParen;}
      bool isOpenBrace() const {return type_==TT_OpenBrace;}
      bool isCloseBrace() const {return type_==TT_CloseBrace;}
      bool isOpenBracket() const {return type_==TT_OpenBracket;}
      bool isCloseBracket() const {return type_==TT_CloseBracket;}
      bool isComma() const {return type_==TT_Comma;}
      bool isDot() const {return type_==TT_Dot;}
      bool isQuotedString() const {return type_==TT_QuotedString;}

      bool isPlus() const {return type_==TT_Plus;}
      bool isMinus() const {return type_==TT_Minus;}
      bool isTimes() const {return type_==TT_Times;}
      bool isDivide() const {return type_==TT_Divide;}

      bool isConstant() const {return type_==TT_Constant;}
      bool isName() const {return type_==TT_Name;}

      bool isEquality() const {return type_==TT_Equal;}
      bool isLessThan() const {return type_==TT_LessThan;}
      bool isGreaterThan() const {return type_==TT_GreaterThan;}
      bool isLessThanOrEqual() const {return type_==TT_LessThanOrEqual;}
      bool isGreaterThanOrEqual() const {return type_==TT_GreaterThanOrEqual;}
      bool isAnd() const {return type_==TT_And;}
      bool isOr() const {return type_==TT_Or;}
      bool isNot() const {return type_==TT_Not;}
      bool isNotEqual() const {return type_==TT_NotEqual;}

      bool isEnd() const {return type_==TT_End;}

      double value() const ;
      const string& name() const ;

      const string& tok() const {return tok_;}
      string stripQuotes() const ;
      void print(ostream& os) const {os << tok_;}
    private:
      string tok_;
      TokenType type_;
    };
}

namespace std
{
  inline ostream& operator<<(ostream& os, const SundanceXML::Token& t)
    {
      t.print(os);
      return os;
    }
}

namespace Teuchos
{

  inline std::string toString(const SundanceXML::Token& tok)
    {
      return tok.tok();
    }

}
#endif
