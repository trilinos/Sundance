#include "SundanceExprParser.hpp"
#include "SundanceExceptions.hpp"


using namespace SundanceUtils;
using namespace SundanceXML;
using namespace Teuchos;


XMLObject ExprParser::evaluateAssignment(const string& line) 
{
	XMLObject result;
  ExprScanner scanner(line);
  assignLevel(scanner, result);
	return result;
}

XMLObject ExprParser::evaluate(const string& line)
{
	XMLObject result;
  ExprScanner scanner(line);
  addLevel(scanner, result);
	return result;
}

void ExprParser::assignLevel(const ExprScanner& scanner, XMLObject& result)
{
	string tmp;
  //  cerr << "entering assign level " << endl;

  if (scanner.peekAtNext().isAssign())
    {
      Token token = scanner.pop();
      tmp = token.tok();
      result = XMLObject("Assign");
      result.addAttribute("target", tmp);
      token = scanner.pop();
      XMLObject lhs;
      addLevel(scanner, lhs);
      result.addChild(lhs);
    }
  else
    {
      addLevel(scanner, result);
    }

  //  cerr << "leaving assign level xml=" << result << endl;
}


void ExprParser::addLevel(const ExprScanner& scanner,
                          XMLObject& result) 
{
  XMLObject tmp;
  bool isPlus;
  bool isFirstTerm = true;
  bool multipleTerms = false;

  //  cerr << "entering add level " << endl;
  multLevel(scanner, tmp);
  Token token = scanner.peek();

  while ((isPlus=token.isPlus()) || token.isMinus())
    {
      token = scanner.pop();
      multipleTerms = true;
      if (isFirstTerm) 
        { 
          result = XMLObject("Sum");
          result.addChild(tmp);
          isFirstTerm = false;
        }						
      multLevel(scanner, tmp);
      token = scanner.peek();					
      if (isPlus) 
        {
          result.addChild(tmp);
        }
      else
        {
          XMLObject term("UnaryMinus");
          term.addChild(tmp);
          result.addChild(term);
        }
    };
  if (isFirstTerm) result = tmp;
  //  cerr << "leaving add level xml=" << result << endl;
}
				
void ExprParser::multLevel(const ExprScanner& scanner,
													 XMLObject& result) 
{
  //  cerr << "entering mult level " << endl;
  XMLObject tmp;
  bool isTimes;
  bool isFirstTerm = true;
  bool multipleTerms = false;

  unaryLevel(scanner, tmp);
  Token token = scanner.peek();

  while ((isTimes=token.isTimes()) || token.isDivide())
    {
      token = scanner.pop();
      multipleTerms = true;
      if (isFirstTerm) 
        { 
          result = XMLObject("Product");
          result.addChild(tmp);
          isFirstTerm = false;
        }						
      unaryLevel(scanner, tmp);
      token = scanner.peek();					
      if (isTimes) 
        {
          result.addChild(tmp);
        }
      else
        {
          XMLObject term("Reciprocal");
          term.addChild(tmp);
          result.addChild(term);
        }
    };
  if (isFirstTerm) result = tmp;
  //  cerr << "leaving mult level xml=" << result << endl;
}
				
void ExprParser::unaryLevel(const ExprScanner& scanner,
                            XMLObject& result)
{
  //  cerr << "entering unary level " << endl;
  XMLObject tmp;
  bool isMinus=false;
  bool isFunction=false;
  bool isIndexed=false;

  // look at the first token, to see if it is a function name, a unary
  // minus, or a unary plus.
  Token start = scanner.peek();
  Token token;

  if (start.isName())
    {
      // if the next token is an open paren, the initial name is to be
      // interpreted as a function.
      if (scanner.peekAtNext().isOpenParen())
        {
          isFunction = true;
          token = scanner.pop();
        }
      // if the next token is an open bracket, the initial name is to be
      // interpreted as a function.
      if (scanner.peekAtNext().isOpenBracket())
        {
          isIndexed = true;
          token = scanner.pop();
        }
    }
  else if (start.isMinus())
    {
      isMinus = true;
      token = scanner.pop();
    }
  else if (start.isPlus())
    {
      token = scanner.pop();
    }
			
  parenLevel(scanner, tmp);
			
  if (isMinus)
    {
      result = XMLObject("UnaryMinus");
      result.addChild(tmp);
    }
  else if (isFunction) 
    {
      result = XMLObject("Function");
      result.addAttribute("name", start.name());
      if (tmp.getTag() != "ArgumentList")
        {
          result.addChild(tmp);
        }
      else
        {
          for (int i=0; i<tmp.numChildren(); i++)
            {
              result.addChild(tmp.getChild(i));
            }
        }
    }
  else if (isIndexed) 
    {
      result = XMLObject("IndexedVar");
      result.addAttribute("name", start.name());
      if (tmp.getTag() != "ArgumentList")
        {
          result.addChild(tmp);
        }
      else
        {
          for (int i=0; i<tmp.numChildren(); i++)
            {
              result.addChild(tmp.getChild(i));
            }
        }
    }
  else result = tmp;

  token = scanner.peek();

  TEST_FOR_EXCEPTION(token.isOpenParen(), RuntimeError,
                     "open paren at unary level: " << scanner.showError());

  // cerr << "leaving unary level xml=" << result << endl;
}

void ExprParser::parenLevel(const ExprScanner& scanner,
                            XMLObject& result)
{
  // cerr << "entering paren level" << endl;
  XMLObject tmp;
			
  Token token = scanner.peek();
			
  if (token.isOpenParen())
    {
      token = scanner.pop();
      result = XMLObject("ArgumentList");
      TEST_FOR_EXCEPTION(token.isCloseParen(), RuntimeError,
                         "empty paren: " << scanner.showError());
      do
        {
          if (token.isComma()) 
            {
              token = scanner.pop();
            }
          addLevel(scanner, tmp);
          token = scanner.peek();
          result.addChild(tmp);
        } 
      while (token.isComma());
      // expect end paren. 
      TEST_FOR_EXCEPTION(!token.isCloseParen(), RuntimeError,
                         "close paren not found: " <<  scanner.showError());
      token = scanner.pop();
      if (result.numChildren()==1) 
        {
          tmp = result.getChild(0);
          result = tmp;
        }
    }
  else if (token.isOpenBrace())
    {
      result = XMLObject("List");
      token = scanner.pop();
      TEST_FOR_EXCEPTION(token.isCloseBrace(), RuntimeError,
                         "empty braces: " << scanner.showError());
      do
        {
          if (token.isComma()) 
            token = scanner.pop();
          addLevel(scanner, tmp);
          token = scanner.peek();
          result.addChild(tmp);
        } 
      while (token.isComma());
      // expect end brace. 
      TEST_FOR_EXCEPTION(!token.isCloseBrace(), RuntimeError,
                         "close brace not found: " <<  scanner.showError());
      token = scanner.pop();
    }
  else if (token.isOpenBracket())
    {
      result = XMLObject("Index");
      token = scanner.pop();
      TEST_FOR_EXCEPTION(token.isCloseBracket(), RuntimeError,
                         "empty bracket: " << scanner.showError());
      do
        {
          if (token.isComma()) 
            token = scanner.pop();
          addLevel(scanner, tmp);
          token = scanner.peek();
          result.addChild(tmp);
        } 
      while (token.isComma());
      // expect end brace. 
      TEST_FOR_EXCEPTION(!token.isCloseBracket(), RuntimeError,
                         "close bracket not found: " <<  scanner.showError());
      token = scanner.pop();
    }
  else 
    {
      primitiveLevel(scanner, result);
      token = scanner.pop();
    }

  //  cerr << "leaving paren level xml=" << result << endl;
}


void ExprParser::primitiveLevel(const ExprScanner& scanner,
                                XMLObject& result) 
{
  //  cerr << "entering paren level" << endl;
  Token token = scanner.peek();
  if (token.isConstant())
    {
      result = XMLObject("Constant");
      result.addAttribute("value", Teuchos::toString(token.value()));
    }
  else if (token.isQuotedString())
    {
      result = XMLObject("String");
      result.addAttribute("value", token.stripQuotes());
    }
  else 
    {
      result = XMLObject("Var");
      result.addAttribute("name", token.name());
    }
  // cerr << "leaving primitive level xml=" << result << endl;
}







