/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SUMEXPR_H
#define SUNDANCE_SUMEXPR_H

#include "SundanceBinaryExpr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using std::string;

  namespace Internal
    {
      /**
       * SumExpr is the internal representation of an addition or subtraction
       * node in the expression tree.
       *
       * Evaluation of derivatives of a sum is given by the sum rule
       * of calculus: the derivative of a sum is the sum of the derivatives.
       * Simple enough. The complication in practice is working with
       * the sparse storage format for the derivatives in an efficient way.
       * Recall from the documentation for EvaluatableExpr and EvalResultSet
       * that entries for derivatives will be created
       * in the results array in the order the requests for derivatives
       * are received during the postprocessing stage, so that the ordering
       * of the results array is essentially arbitrary.
       * Furthermore,
       * a nonzero derivative of the sum might nonetheless
       * be structurally zero in
       * one of the two operands, so it's possible
       * that some derivatives will appear
       * in only one of the two results arrays.
       *\note{There is an ordering for derivatives, of course, but it is of little use here because the left and right sets of nonzero derivatives can be different.}
       *
       * The working out of which left results need to be added to which
       * right results is easily, but slowly, accomplished through
       * derivative-keyed map lookups in the result sets. However, this
       * process can also be done once and for all up front, letting us
       * build a set of arrays giving a set of instructions of
       * which entries should be added. We can store "sum rules"
       * as triplets of integers: one giving the index of the "target"
       * result entry, and two giving the indices of the "source" results
       * in the left and right result tables. If one of the two sources
       * is structurally zero, we can mark that entry with an impossible
       * index; we use an index -1 to mark zero entries.
       *
       * As an example, suppose we are doing the sum \f$F=L+R\f$, and
       * the left and right leaves have the following
       * nonzero derivatives, listed in the left and right results
       * arrays in the order shown:
       * <center>
       * <h4> Left result array </h4>
       * <table>
       * <tr>
       * <td> 0 </td> <td> 1 </td><td> 2 </td> <td> 3 </td><td> 4 </td>
       * </tr>
       * <tr><td> \f$L_u\f$ </td>
       * <td>\f$ L_x\f$ </td>
       * <td>\f$ L\f$ </td>
       * <td>\f$ L_{u_x}\f$ </td>
       * <td>\f$ L_v\f$ </td>
       * </tr></table>
       * </center>
       * <center>
       * <h4> Right result array </h4>
       * <table>
       * <tr>
       * <td> 0 </td> <td> 1 </td><td> 2 </td><td> 3 </td>
       * </tr>
       * <tr><td> \f$R_v\f$ </td>
       * <td>\f$ R\f$ </td>
       * <td>\f$ R_u\f$ </td>
       * <td>\f$ R_{v_x}\f$ </td>
       * </tr></table>
       * </center>
       *
       * Clearly the set of nonzero derivatives of \f$F\f$ is the
       * union of sets of nonzero derivatives of \f$L\f$ and \f$R\f$,
       * or {\f$F, F_u, F_v, F_{u_x}, F_x, F_{v_x}\f$}.
       * For each of these, we can find the pairs of source indices, and
       * build a table of
       * <tt>{target index, left source index, right source index}</tt>
       * triplets.
       * <center>
       * <h4> Sum rule table </h4>
       * <table> <tr>
       * <td> Result index </td> <td> Derivative computed </td>
       * <td> Index into left results </td>
       * <td> Index into right results </td>
       * </tr>
       * <tr><td> 0 </td><td> \f$F \f$ </td><td> 2 </td><td> 1 </td></tr>
       * <tr><td> 1 </td><td> \f$F_u \f$ </td><td> 0 </td><td> 2 </td></tr>
       * <tr><td> 2 </td><td> \f$F_v \f$ </td><td> 4 </td><td> 0 </td></tr>
       * <tr><td> 3 </td><td> \f$F_{u_x} \f$ </td><td> 3 </td><td> -1</td></tr>
       * <tr><td> 4 </td><td> \f$F_{x} \f$ </td><td> 1 </td><td> -1 </td></tr>
       * <tr><td> 5 </td><td> \f$F_{v_x} \f$ </td><td> -1 </td><td> 3</td></tr>
       * </table>
       * </center>
       <br>
       <br>
       * This is then processed with the following operations given in
       * pseudocode:
       * <pre>
       * result[0] = leftResult[2] + rightResult[1];
       *
       * result[1] = leftResult[0] + rightResult[2];
       *
       * result[2] = leftResult[4] + rightResult[0];
       *
       * result[3] = leftResult[3];
       *
       * result[4] = leftResult[1];
       *
       * result[5] =                 rightResult[3];
       * </pre>
       */
      class SumExpr : public BinaryExpr
        {
        public:
          /** */
          SumExpr(const RefCountPtr<ScalarExpr>& a, 
                  const RefCountPtr<ScalarExpr>& b, int sign);

          /** */
          virtual ~SumExpr() {;}

          /** Indicate whether this expression is a "hungry"
           * differential operator that is awaiting an argument. */
          virtual bool isHungryDiffOp() const ;

          /**
           * Indicate whether the given functional derivative is nonzero
           */
          virtual bool hasNonzeroDeriv(const MultipleDeriv& f) const ;


          /** Return the set of derivatives required by the operands
           * of this expression given that this expression
           * requires the set d. For all expressions other than DiffOp,
           * the operand derivative set is identical to the input derivative
           * set. DiffOp will require a different set of derivatives from
           * its operand. */
          virtual Array<DerivSet>
          derivsRequiredFromOperands(const DerivSet& d) const ;

          /** Return true if both left and right return true */
          virtual bool allTermsHaveTestFunctions() const ;


          /** */
          virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

        protected:
          /** */
          virtual bool parenthesizeSelf() const {return true;}
          /** */
          virtual bool parenthesizeOperands() const {return false;}
          /** */
          virtual const string& xmlTag() const ;
          /** */
          virtual const string& opChar() const ;

        private:


        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
