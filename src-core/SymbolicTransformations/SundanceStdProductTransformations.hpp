/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_STDPRODUCTTRANSFORMATION_H
#define SUNDANCE_STDPRODUCTTRANSFORMATION_H

#include "SundanceDefs.hpp"
#include "SundanceProductTransformationSequence.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY 

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore::Internal;

  using std::string;
  using std::ostream;

  namespace Internal
  {
    /**
     * Apply a standard set of transformations
     */
    class StdProductTransformations : public ProductTransformationSequence
    {
    public:
      StdProductTransformations();

      virtual ~StdProductTransformations(){;}
    };
    
    /** 
     * Transform a product by removing a zero term: 
     * \f[
     * x \times 0 \rightarrow 0. 
     * \f]
     * \f[
     * 0 \times x \rightarrow 0. 
     * \f]
     */
    class RemoveZeroFromProduct : public ProductTransformation
    {
    public:
      /** */
      RemoveZeroFromProduct() : ProductTransformation() {;}

      /** */
      virtual ~RemoveZeroFromProduct(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };

    /** 
     * Multiply two constant exprs without transformation 
     */
    class MultiplyConstants : public ProductTransformation
    {
    public:
      /** */
      MultiplyConstants() : ProductTransformation() {;}

      /** */
      virtual ~MultiplyConstants(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };

    /** 
     * Transform a product by moving any constants to the left:
     * \f[
     * x \times a \rightarrow a \times x
     * \f]
     **/
    class MoveConstantsToLeftOfProduct : public ProductTransformation
    {
    public:
      /** */
      MoveConstantsToLeftOfProduct() : ProductTransformation() {;}

      /** */
      virtual ~MoveConstantsToLeftOfProduct(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };

    /**
     * Transform a product by associating any hungry diff op in the
     * left operand with the right operand:
     * \f[
     * (u D_x) v \rightarrow u D_x u
     * \f]
     */
    class AssociateHungryDiffOpWithOperand : public ProductTransformation
    {
    public:
      /** */
      AssociateHungryDiffOpWithOperand() : ProductTransformation() {;}

      /** */
      virtual ~AssociateHungryDiffOpWithOperand(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };

    /**
     * Kill a diff op acting on a constant
     * \f[
     * D_x \alpha \rightarrow 0
     * \f]
     * \f[
     * D_x (\alpha + u) \rightarrow D_x u
     * \f]
     */
    class KillDiffOpOnConstant : public ProductTransformation
    {
    public:
      /** */
      KillDiffOpOnConstant() : ProductTransformation() {;}

      /** */
      virtual ~KillDiffOpOnConstant(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };

    /**
     * Bring a constant outside a diff op
     * \f[
     * D_x (\alpha u) \rightarrow \alpha D_x u
     * \f]
     */
    class BringConstantOutsideDiffOp : public ProductTransformation
    {
    public:
      /** */
      BringConstantOutsideDiffOp() : ProductTransformation() {;}

      /** */
      virtual ~BringConstantOutsideDiffOp(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };
    
    /**
     * Distribute a sum of diff ops over their operand
     * \f[
     * (D_1 + D_2) u \rightarrow D_1 u + D_2 u
     * \f]
     */
    class DistributeSumOfDiffOps : public ProductTransformation
    {
    public:
      /** */
      DistributeSumOfDiffOps() : ProductTransformation() {;}

      /** */
      virtual ~DistributeSumOfDiffOps(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };

    /**
     * Apply a simple diff op
     */
    class ApplySimpleDiffOp : public ProductTransformation
    {
    public:
      /** */
      ApplySimpleDiffOp() : ProductTransformation() {;}

      /** */
      virtual ~ApplySimpleDiffOp(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };

    /** 
     * Rearrange a product whose right operand is 
     * a product including a constant
     * such that constants are grouped on the left:
     * \f[
     * \alpha (\beta u) \rightarrow (\alpha\beta) u
     * \f]
     * \f[
     * u (\alpha v) \rightarrow \alpha (u v)
     * \f]
     **/
    class RearrangeRightProductWithConstant : public ProductTransformation
    {
    public:
      /** */
      RearrangeRightProductWithConstant() : ProductTransformation() {;}

      /** */
      virtual ~RearrangeRightProductWithConstant(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };

    /** 
     * Rearrange a product whose left operand is a product including a constant
     * such that constants are grouped on the left:
     * \f[
     * (\alpha u)\beta \rightarrow \alpha\beta u
     * \f]
     * \f[
     * (\alpha u)v \rightarrow \alpha (u v)
     * \f]
     **/
    class RearrangeLeftProductWithConstant : public ProductTransformation
    {
    public:
      /** */
      RearrangeLeftProductWithConstant() : ProductTransformation() {;}

      /** */
      virtual ~RearrangeLeftProductWithConstant(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };

    /** 
     * Rearrange a product of a constant and an integral so that
     * the constant is under the integral sign:
     * \f[
     * \alpha \int u \rightarrow \int \alpha u
     * \f]
     **/
    class TakeConstantUnderIntegralSign : public ProductTransformation
    {
    public:
      /** */
      TakeConstantUnderIntegralSign() : ProductTransformation() {;}

      /** */
      virtual ~TakeConstantUnderIntegralSign(){;}

      /** */
      virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                               RefCountPtr<ScalarExpr>& rtn) const ;
    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif