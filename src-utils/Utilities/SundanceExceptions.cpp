/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExceptions.hpp"

using namespace SundanceUtils;

InternalError::InternalError(const string& msg)
  : std::logic_error(msg)
{;}

RuntimeError::RuntimeError(const string& msg)
  : std::runtime_error(msg)
{;}

BadSymbolicsError::BadSymbolicsError(const string& msg)
  : RuntimeError(msg)
{;}


