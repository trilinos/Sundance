/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExpr.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceStdSumTransformations.hpp"

#include "SundanceSumExpr.hpp"
#include "SundanceProductExpr.hpp"
#include "SundanceConstantExpr.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceDiffOp.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceUnaryMinus.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSumOfBCs.hpp"
#include "SundanceNullCellFilterStub.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace SundanceCore::Internal;


StdSumTransformations::StdSumTransformations()
  : SumTransformationSequence()
{
  append(rcp(new RemoveZeroFromSum()));
  append(rcp(new SumConstants()));
  append(rcp(new MoveConstantsToLeftOfSum()));
  append(rcp(new RearrangeRightSumWithConstant()));
  append(rcp(new RearrangeLeftSumWithConstant()));
  append(rcp(new SumIntegrals()));
}

bool RemoveZeroFromSum::doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                                    int sign, RefCountPtr<ScalarExpr>& rtn) const
{
  /* Check for the trivial case of operation with zero */
  
  /* If I'm constant and my value is zero, return other */
  const ConstantExpr* cl = dynamic_cast<const ConstantExpr*>(left.get());
  if (cl != 0 && cl->value()==0.0)
    {
      if (verbosity() > 1)
        {
          Out::println("RemoveZeroFromSum identified left operand as zero.");
          Out::println("Applying transformation 0 + x --> x");
        }
      rtn = chooseSign(sign, right);
      return true;
    }

  /* If other is zero, return me */
  const ConstantExpr* cr = dynamic_cast<const ConstantExpr*>(right.get());
  if (cr != 0 && cr->value()==0.0) 
    {
      if (verbosity() > 1)
        {
          Out::println("RemoveZeroFromSum identified right operand as zero.");
          Out::println("Applying transformation x + 0 --> x");
        }
      rtn = left;
      return true;
    }

  /* otherwise, do no transformation */
  return false;
}

bool MoveConstantsToLeftOfSum::doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                                      int sign, RefCountPtr<ScalarExpr>& rtn) const
{
  /* if the right operand is a constant, 
   * transform u +/- alpha --> +/- alpha + u */
  if (right->isConstant())
    {
      if (verbosity() > 1)
        {
          Out::println("MoveConstantsToLeftOfSum identified right "
                       "operand as constant.");
        }
      rtn = getScalar(Expr::handle(chooseSign(sign, right)) 
                      + Expr::handle(left));
      return true;
    }

  return false;
}

bool SumConstants::doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               int sign, RefCountPtr<ScalarExpr>& rtn) const
{
  /* Check to see if both are constant. If so, sum them and return */
  if (left->isConstant() && right->isConstant())
    {
      if (verbosity() > 1)
        {
          Out::println("SumConstants identified both "
                       "operands as constant. No transformations applied.");
        }
      rtn = rcp(new SumExpr(left, right, sign));
      return true;
    }
  return false;
}

bool RearrangeRightSumWithConstant::doTransform(const RefCountPtr<ScalarExpr>& left, 
                                                const RefCountPtr<ScalarExpr>& right,
                                                int sign, RefCountPtr<ScalarExpr>& rtn) const
{
  const SumExpr* sRight = dynamic_cast<const SumExpr*>(right.get());

  if (sRight != 0)
    {
      /* The case in which the right operand of right is a constant
       * should have been transformed away by now. Do a paranoid 
       * check to make sure this hasn't happened */
      TEST_FOR_EXCEPTION(sRight->rightScalar()->isConstant(),
                         InternalError,
                         "RearrangeRightSumWithConstant: unexpected case, "
                         "constant expr"
                         << sRight->right() << " found as right operand "
                         "in sum " << right->toString());

      if (sRight->leftScalar()->isConstant())
        {      
          /* If left operand is a constant, transform
           * alpha + s1*(beta + s2*u) --> (alpha + s1*beta) + s1*s2*u */
          if (left->isConstant())
            {
              if (verbosity() > 1)
                {
                  Out::println("RearrangeRightSumWithConstant::doTransform "
                               "identified right "
                               "operand as sum involving a constant, "
                               "and left operand as a constant. Applying "
                               "transformation alpha + (beta+u) "
                               "--> (alpha + beta) + u.");
                }
              int s1 = sign;
              int s2 = sRight->sign();
              Expr alpha = Expr::handle(left);
              Expr beta = sRight->left();
              Expr u = sRight->right();
              rtn = getScalar((alpha + chooseSign(s1, beta)) + chooseSign(s1*s2,u));
            }
          else  /* if left operand is non-constant, transform
                 * u + s1*(alpha + s2*v) --> s1*alpha + (u + s1*s2*v) */
            {
              if (verbosity() > 1)
                {
                  Out::println("RearrangeRightSumWithConstant::doTransform "
                               "identified right "
                               "operand as sum involving a constant, "
                               "and left operand as non-constant. Applying "
                               "transformation u + (alpha + v) "
                               "--> alpha + (u + v)");
                }
              int s1 = sign;
              int s2 = sRight->sign();
              Expr u = Expr::handle(left);
              Expr alpha = sRight->left();
              Expr v = sRight->right();
              rtn = getScalar(chooseSign(s1, alpha) + (u + chooseSign(s1*s2, v)));
            }
          return true;
        }
    }
  return false;
}


bool RearrangeLeftSumWithConstant::doTransform(const RefCountPtr<ScalarExpr>& left, 
                                                const RefCountPtr<ScalarExpr>& right,
                                                int sign, RefCountPtr<ScalarExpr>& rtn) const
{
  const SumExpr* sLeft = dynamic_cast<const SumExpr*>(left.get());

  if (sLeft != 0 && !left->isConstant())
    {
      /* The case in which the right operand of left is a constant
       * should have been transformed away by now. Do a paranoid 
       * check to make sure this hasn't happened */
      TEST_FOR_EXCEPTION(sLeft->rightScalar()->isConstant(),
                         InternalError,
                         "RearrangeLeftSumWithConstant::doTransform "
                         ": unexpected case, constant expr"
                         << sLeft->right() << " found as right operand "
                         "in sum " << left->toString());
      
      if (sLeft->leftScalar()->isConstant())
        {
          /* If right operand is a constant, transform 
           * (alpha + s1*u) + s2*beta --> (alpha + s2*beta) + s1*u */
          if (right->isConstant())
            {
              if (verbosity() > 1)
                {
                  Out::println("RearrangeLeftSumWithConstant::doTransform "
                               "identified right "
                               "operand as constant, "
                               "and left operand as sum involving "
                               "a constant. Applying "
                               "transformation (alpha + u) + beta "
                               "--> (alpha + beta) + u.");
                }
              int s2 = sign;
              int s1 = sLeft->sign();
              Expr alpha = sLeft->left();
              Expr beta = Expr::handle(right);
              Expr u = sLeft->right();
              rtn =  getScalar((alpha + chooseSign(s2, beta)) + chooseSign(s1, u));
            }
          else /* if right operand is a non-constant, transform 
                * (alpha + s1*u) + s2*v --> alpha + (s1*u + s2*v) */
            {
              if (verbosity() > 1)
                {
                  Out::println("RearrangeLeftSumWithConstant::doTransform "
                               "identified right "
                               "operand as non-constant, "
                               "and left operand as sum involving "
                               "a constant. Applying "
                               "transformation (alpha + u) + v "
                               "--> alpha + (u + v).");
                }
              int s2 = sign;
              int s1 = sLeft->sign();
              Expr alpha = sLeft->left();
              Expr u = sLeft->right();
              Expr v = Expr::handle(right);
              rtn =  getScalar(alpha 
                + (chooseSign(s1, u) + chooseSign(s2, v)));
            }
          return true;
        }
    }
  return false;
}

bool SumIntegrals::doTransform(const RefCountPtr<ScalarExpr>& left, 
                               const RefCountPtr<ScalarExpr>& right,
                               int sign, RefCountPtr<ScalarExpr>& rtn) const
{
  const SumOfIntegrals* sLeft 
    = dynamic_cast<const SumOfIntegrals*>(left.get());
  const SumOfIntegrals* sRight 
    = dynamic_cast<const SumOfIntegrals*>(right.get());

  if (sLeft != 0 || sRight != 0)
    {
      /* make sure we don't have a case where one is an essential BC and
       * the other is not */
      bool leftIsBC = (dynamic_cast<const SumOfBCs*>(sLeft) != 0);
      bool rightIsBC = (dynamic_cast<const SumOfBCs*>(sRight) != 0);
      TEST_FOR_EXCEPTION((leftIsBC && !rightIsBC)
                         || (!leftIsBC && rightIsBC), RuntimeError,
                         "Attempting to add EssentialBC and non-EssentialBC "
                         "integrals: L=" << left->toString() << ", R="
                         << right->toString());

      if (sLeft != 0 && sRight != 0)
        {
          SumOfIntegrals* l;
          if (!leftIsBC) l = new SumOfIntegrals(*sLeft);
          else l = new SumOfBCs(*dynamic_cast<const SumOfBCs*>(sLeft));
          l->merge(sRight, sign);
          rtn = rcp(l);
          return true;
        }

      /* at this point, one of the terms is a global equation. BCs should
       * not be involved at this point */
      TEST_FOR_EXCEPTION(leftIsBC, RuntimeError,
                         "Attempting to add a BC " << left->toString()
                         << " and a global expression " << right->toString());

      TEST_FOR_EXCEPTION(rightIsBC, RuntimeError,
                         "Attempting to add a BC " << right->toString()
                         << " and a global expression " << left->toString());

      if (sLeft != 0 && sRight == 0)
        {
          SumOfIntegrals* l = new SumOfIntegrals(*sLeft);
          const SpatiallyConstantExpr* cRight 
            = dynamic_cast<const SpatiallyConstantExpr*>(right.get());

          TEST_FOR_EXCEPTION(cRight == 0, InternalError,
                             "Attempting to add non-constant expression "
                             << right->toString() << " to an integral");

          Expr r = Integral(l->nullRegion(), Expr::handle(right));
          const SumOfIntegrals* sr 
            = dynamic_cast<const SumOfIntegrals*>(r.ptr().get());
          l->merge(sr, sign);
          rtn = rcp(l);
          return true;
        }
      else
        {
          SumOfIntegrals* r = new SumOfIntegrals(*sRight);
          if (sign < 0) r->changeSign();

          const SpatiallyConstantExpr* cLeft 
            = dynamic_cast<const SpatiallyConstantExpr*>(left.get());

          TEST_FOR_EXCEPTION(cLeft == 0, InternalError,
                             "Attempting to add non-constant expression "
                             << left->toString() << " to an integral");

          Expr l = Integral(r->nullRegion(), Expr::handle(right));
          const SumOfIntegrals* sl 
            = dynamic_cast<const SumOfIntegrals*>(l.ptr().get());
          r->merge(sl, 1);
          rtn = rcp(r);
          return true;
        }

    }
  else
    {
      return false;
    }
  
}
