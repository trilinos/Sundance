/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFunctorDomain.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceUtils;
using namespace Teuchos;


FunctorDomain::FunctorDomain() {;}

double FunctorDomain::lowerBound() const
{
  TEST_FOR_EXCEPTION(true, InternalError, 
                     "FunctorDomain::lowerBound() called for a domain without "
                     "a lower bound");
  return 0.0;
}

double FunctorDomain::upperBound() const
{
  TEST_FOR_EXCEPTION(true, InternalError, 
                     "FunctorDomain::upperBound() called for a domain without "
                     "an upper bound");

  return 0.0;
}

double FunctorDomain::excludedPoint() const
{
  TEST_FOR_EXCEPTION(true, InternalError, 
                     "FunctorDomain::excludedPoint() called for a domain without "
                     "an excluded point");

  return 0.0;
}

UnboundedDomain::UnboundedDomain() {;}

PositiveDomain::PositiveDomain() {;}

BoundedDomain::BoundedDomain(const double& lower, const double& upper) 
  : FunctorDomain(), lower_(lower), upper_(upper) 
{;}

LowerBoundedDomain::LowerBoundedDomain(const double& lower)
  : FunctorDomain(), lower_(lower)
{;}

NonzeroDomain::NonzeroDomain() {;}

