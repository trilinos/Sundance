#include "SundanceExprScanner.hpp"
#include "SundanceExceptions.hpp"


using namespace SundanceXML;
using namespace SundanceUtils;
using namespace Teuchos;


ExprScanner::ExprScanner(const string& line)
	: tok_(0), ptr_(0)
{
	Token tok;
	int place = 0;

	while (!(tok=getToken(line, place)).isEnd())
    {
      tok_.append(tok);
    }
	tok_.append(tok);
}

bool ExprScanner::hasMore() const 
{
	return (ptr_ < tok_.length());
}

const Token& ExprScanner::pop() const
{
  const Token& rtn = tok_[ptr_];
  ptr_++;
  return rtn;

	return tok_[0];
}

const Token& ExprScanner::peek() const
{
  const Token& rtn = tok_[ptr_];
  return rtn;
}

const Token& ExprScanner::peekAtNext() const
{
  const Token& rtn = tok_[ptr_+1];
  return rtn;
}

string ExprScanner::showError() const 
{
	string rtn;
	for (int p=0; p<ptr_; p++)
		{
			rtn += tok_[p].tok();
		}
	rtn += "\n<ERROR>\n";
	for (int p=ptr_; p<tok_.length(); p++)
		{
			rtn += tok_[p].tok();
		}
  rtn += "\n</ERROR>\n";
	return rtn;
}
 
bool ExprScanner::isWhite(char c) const 
{
  if (c==' ' || c=='\t' || c=='\n' || c=='\r') return true;
  return false;
}

bool ExprScanner::isDelimiter(char c) const 
{
  if (strchr("()+-*/,{};", c) && c!='\0') return true;
  return false;
}

bool ExprScanner::isLogical(char c) const 
{
  if (strchr("=!<>&|", c) && c!='\0') return true;
  return false;
}

Token ExprScanner::getToken(const string& line, int& place) const 
{
  string currentToken;
			
  // skip whitespace
  while (place < (int) line.length() && isWhite(line[place])) 
    {
      place++;
    }
			
  // check for end of string
  if (place >= (int) line.length() || line[place]=='\0') 
    {
      return Token(";");
    }
			
			
  // check for delimiter
  if (isDelimiter(line[place]))
    {
      currentToken = line[place];
      place++;
      return Token(currentToken);
    }

  // check for logical operator. Some logicals are two chars long. 
  if (isLogical(line[place]))
    {
      currentToken = line[place];
      if (isLogical(line[place+1]))
        {
          currentToken += line[place+1];
          place++;
        }
      place++;
      return Token(currentToken);
    }
			
  // check for numerical value
			
  if (isdigit(line[place]))
    {
      currentToken = line[place];
      place++;
      while (!(place >= (int) line.length() 
               || isDelimiter(line[place]) 
               || isWhite(line[place])))
        {
          currentToken += line[place];
          place++;
        }
      return Token(currentToken);
    }
			
  // check for named variable
  if (isalpha(line[place]))
    {
      currentToken = line[place];
      place++;
      while (!(place >= (int) line.length() 
               || isDelimiter(line[place]) || isWhite(line[place])))
        {
          currentToken += line[place];
          place++;
        }
      return Token(currentToken);
    }
			
  // haven't found anything sensible...
  TEST_FOR_EXCEPTION(true, RuntimeError,
                     "undecipherable line " << line);

	return Token(""); // -Wall
}


			
			
	

