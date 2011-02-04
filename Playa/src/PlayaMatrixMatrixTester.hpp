/* @HEADER@ */
//   
 /* @HEADER@ */


#ifndef PLAYA_MATRIXMATRIXTESTER_HPP
#define PLAYA_MATRIXMATRIXTESTER_HPP

#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaEpetraMatrixMatrixProduct.hpp"
#include "PlayaEpetraMatrixMatrixSum.hpp"
#include "PlayaEpetraMatrixOps.hpp"
#include "PlayaEpetraMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "PlayaLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaSimpleComposedOpImpl.hpp"
#include "PlayaSimpleScaledOpImpl.hpp"
#include "PlayaSimpleAddedOpImpl.hpp"
#include "PlayaSimpleDiagonalOpImpl.hpp"
#include "PlayaRandomSparseMatrixBuilderImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;



namespace Playa
{

/** */
template <class Scalar>
class MatrixMatrixTester : public TesterBase<Scalar>
{
public:
  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  MatrixMatrixTester(const LinearOperator<Scalar>& A,
    const LinearOperator<Scalar>& B,
    const TestSpecifier<Scalar>& sumSpec,
    const TestSpecifier<Scalar>& prodSpec,
    const TestSpecifier<Scalar>& diagSpec,
    const TestSpecifier<Scalar>& diagLeftProdSpec,
    const TestSpecifier<Scalar>& diagRightProdSpec);

  /** */
  bool runAllTests() const ;

  /** */
  bool sumTest() const ;

  /** */
  bool prodTest() const ;

  /** */
  bool diagTest() const ;

  /** */
  bool diagLeftProdTest() const ;

  /** */
  bool diagRightProdTest() const ;


private:

  LinearOperator<Scalar> A_;

  LinearOperator<Scalar> B_;

  TestSpecifier<Scalar> sumSpec_;

  TestSpecifier<Scalar> prodSpec_;

  TestSpecifier<Scalar> diagSpec_;

  TestSpecifier<Scalar> diagLeftProdSpec_;

  TestSpecifier<Scalar> diagRightProdSpec_;

};

template <class Scalar> 
inline MatrixMatrixTester<Scalar>
::MatrixMatrixTester(const LinearOperator<Scalar>& A,
  const LinearOperator<Scalar>& B,
  const TestSpecifier<Scalar>& sumSpec,
  const TestSpecifier<Scalar>& prodSpec,
  const TestSpecifier<Scalar>& diagSpec,
  const TestSpecifier<Scalar>& diagRightProdSpec,
  const TestSpecifier<Scalar>& diagLeftProdSpec)
  : TesterBase<Scalar>(), 
    A_(A),
    B_(B),
    sumSpec_(sumSpec),
    prodSpec_(prodSpec),
    diagSpec_(diagSpec),
    diagLeftProdSpec_(diagLeftProdSpec),
    diagRightProdSpec_(diagRightProdSpec)
{;}

template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::runAllTests() const
{
  bool pass = true;

  pass = this->sumTest() && pass;
  pass = this->prodTest() && pass;
  pass = this->diagTest() && pass;
  pass = this->diagLeftProdTest() && pass;
  pass = this->diagRightProdTest() && pass;

  return pass;
}

template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::diagTest() const 
{
  if (diagSpec_.doTest())
  {
    Vector<Scalar> x = B_.domain().createMember();
    Vector<Scalar> d = B_.domain().createMember();
    randomizeVec(x);
    randomizeVec(d);
    LinearOperator<Scalar> D0 = diagonalOperator(d);
    LinearOperator<Scalar> D = makeEpetraDiagonalMatrix(d);
    Vector<Scalar> d1 = getEpetraDiagonal(D);

    Out::root() << "computing implicit product y1 = D*x..." << std::endl;
    Vector<Scalar> y1 = D0*x;
    Out::root() << "computing explicit product y2 = D*x..." << std::endl;
    Vector<Scalar> y2 = D*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;
    
    Out::root() << "comparing recovered and original diagonals" << std::endl;
    ScalarMag err2 = (d - d1).norm2();
    Out::root() << "|d1-d2| = " << err2 << std::endl;
    
    return checkTest(prodSpec_, err+err2, "matrix-matrix multiply");
    
  }
  Out::root() << "skipping matrix-matrix multiply test..." << std::endl;
  return true;
}



template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::sumTest() const 
{
  if (sumSpec_.doTest())
  {
    /* skip incompatible matrices. This will occur when we're testing
     * multiplication of rectangular matrices */
    if (A_.range() != B_.range() || A_.domain() != B_.domain())
    {
      Out::root() << "skipping sum on incompatible matrices" << std::endl;
      return true;
    }
    /* If here, the sum should work */
    Out::root() << "running matrix-matrix multiply test..." << std::endl;
    LinearOperator<Scalar> implicitAdd = A_ + B_;
    LinearOperator<Scalar> explicitAdd = epetraMatrixMatrixSum(A_, B_);

    Vector<Scalar> x = B_.domain().createMember();
    randomizeVec(x);
    Out::root() << "computing implicit sum y1 = (A+B)*x..." << std::endl;
    Vector<Scalar> y1 = implicitAdd*x;
    Out::root() << "computing explicit sum y2 = (A+B)*x..." << std::endl;
    Vector<Scalar> y2 = explicitAdd*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;
    return checkTest(prodSpec_, err, "matrix-matrix multiply");
    
  }
  Out::root() << "skipping matrix-matrix multiply test..." << std::endl;
  return true;
}


template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::prodTest() const 
{
  if (prodSpec_.doTest())
  {
    Out::root() << "running matrix-matrix multiply test..." << std::endl;
    LinearOperator<Scalar> composed = A_ * B_;
    LinearOperator<Scalar> multiplied = epetraMatrixMatrixProduct(A_, B_);

    Vector<Scalar> x = B_.domain().createMember();
    randomizeVec(x);
    Out::root() << "computing implicit product y1 = (A*B)*x..." << std::endl;
    Vector<Scalar> y1 = composed*x;
    Out::root() << "computing explicit product y2 = (A*B)*x..." << std::endl;
    Vector<Scalar> y2 = multiplied*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;
    return checkTest(prodSpec_, err, "matrix-matrix multiply");
  }
  Out::root() << "skipping matrix-matrix multiply test..." << std::endl;
  return true;
}


template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::diagLeftProdTest() const 
{
  if (diagLeftProdSpec_.doTest())
  {
    Out::root() << "running diagonal*matrix multiplication test..." << std::endl;

    Vector<Scalar> x = A_.domain().createMember();
    randomizeVec(x);

    Vector<Scalar> d = A_.range().createMember();
    randomizeVec(d);
        
    LinearOperator<Scalar> D = diagonalOperator(d);
    LinearOperator<Scalar> DA = epetraLeftScale(d, A_);

    Out::root() << "computing implicit y1 = D*A*x..." << std::endl;
    Vector<Scalar> y1 = D*A_*x;
    Out::root() << "computing explicit y2 = D*A*x..." << std::endl;
    Vector<Scalar> y2 = DA*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;

    return checkTest(diagLeftProdSpec_, err, "diagonal*matrix multiplication");
  }
  Out::root() << "skipping diagonal matrix-matrix test..." << std::endl;
  return true;
}

  
template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::diagRightProdTest() const 
{
  if (diagRightProdSpec_.doTest())
  {
    Out::root() << "running diagonal*matrix multiplication test..." << std::endl;

    Vector<Scalar> x = A_.domain().createMember();
    randomizeVec(x);

    Vector<Scalar> d = A_.domain().createMember();
    randomizeVec(d);
        
    LinearOperator<Scalar> D = diagonalOperator(d);
    LinearOperator<Scalar> AD = epetraRightScale(A_, d);

    Out::root() << "computing implicit y1 = A*D*x..." << std::endl;
    Vector<Scalar> y1 = A_*D*x;
    Out::root() << "computing explicit y2 = A*D*x..." << std::endl;
    Vector<Scalar> y2 = AD*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;

    return checkTest(diagLeftProdSpec_, err, "matrix*diagonal multiplication");
  }
  Out::root() << "skipping diagonal matrix-matrix test..." << std::endl;
  return true;
}

  
  
}
#endif
