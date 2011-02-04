/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaExceptions.hpp"

using namespace Playa;

InternalError::InternalError(const std::string& msg)
  : std::logic_error(msg)
{;}

RuntimeError::RuntimeError(const std::string& msg)
  : std::runtime_error(msg)
{;}

