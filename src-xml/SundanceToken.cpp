#include "SundanceToken.hpp"
#include "SundanceExceptions.hpp"
#include <ctype.h>


using namespace SundanceXML;
using namespace SundanceUtils;
using namespace Teuchos;

Token::Token(const string& tok)
	: tok_(tok),
		type_(TT_Unknown)
{
  TEST_FOR_EXCEPTION(tok.length() <= 0, InternalError,
                     "empty string in Token ctor");

	if (isdigit(tok_[0])) type_=TT_Constant;
	
	else if (isalpha(tok_[0])) type_=TT_Name;
	
	else if (tok_ == "=") type_=TT_Assign;
	else if (tok_[0]=='(') type_=TT_OpenParen;
	else if (tok_[0]==')') type_=TT_CloseParen;
	else if (tok_[0]=='{') type_=TT_OpenBrace;
	else if (tok_[0]=='}') type_=TT_CloseBrace;
	else if (tok_[0]=='[') type_=TT_OpenBracket;
	else if (tok_[0]==']') type_=TT_CloseBracket;
	else if (tok_[0]==',') type_=TT_Comma;
	else if (tok_[0]=='.') type_=TT_Dot;
	else if (tok_[0]=='\'') type_=TT_QuotedString;

	else if (tok_[0]=='+') type_=TT_Plus;
	else if (tok_[0]=='-') type_=TT_Minus;
	else if (tok_[0]=='*') type_=TT_Times;
	else if (tok_[0]=='/') type_=TT_Divide;

	else if (tok_ == "==") type_=TT_Equal;
	else if (tok_ == "!=") type_=TT_NotEqual;
	else if (tok_ == "<=") type_=TT_LessThanOrEqual;
	else if (tok_ == ">=") type_=TT_GreaterThanOrEqual;
	else if (tok_ == ">") type_=TT_GreaterThan;
	else if (tok_ == "<") type_=TT_LessThan;
	else if (tok_ == "!") type_=TT_Not;
	else if (tok_ == "&&") type_=TT_And;
	else if (tok_ == "||") type_=TT_Or;
	else type_=TT_End;
}

double Token::value() const 
{
	double val=0.0;

  TEST_FOR_EXCEPTION(type_!=TT_Constant, InternalError,
                     "Token::value() called for non-Constant type token");

  return atof(tok_.c_str());
}

const string& Token::name() const 
{
  TEST_FOR_EXCEPTION(type_!=TT_Name, InternalError,
                     "Token::name() called for non-Name type token");
	return tok_;
}


string Token::stripQuotes() const 
{
  TEST_FOR_EXCEPTION(type_!=TT_QuotedString, InternalError,
                     "Token::name() called for non-Name type token");
  int n = tok_.length();
  string rtn;
  if (n==2) return rtn;
  return tok_.substr(1, n-2);

}


