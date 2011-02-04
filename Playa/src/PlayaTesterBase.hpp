/* @HEADER@ */
//   
/* @HEADER@ */


#ifndef PLAYA_TESTERBASE_HPP
#define PLAYA_TESTERBASE_HPP

#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaTestSpecifier.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "PlayaOut.hpp"


using namespace Teuchos;



namespace Playa
{

/** */
template <class Scalar>
class TesterBase
{
public:
  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  TesterBase(){;}

  /** */
  virtual ~TesterBase(){;}

  /** */
  virtual bool runAllTests() const = 0 ;


  /** */
  bool checkTest(const TestSpecifier<Scalar>& spec,
    const ScalarMag& err, 
    const std::string& testName) const ;

  /** */
  void randomizeVec(Vector<Scalar>& x) const ;

};

template <class Scalar> 
inline void TesterBase<Scalar>
::randomizeVec(Vector<Scalar>& x) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  x.randomize();
    
}

template <class Scalar> 
inline bool TesterBase<Scalar>
::checkTest(const TestSpecifier<Scalar>& spec,
  const ScalarMag& err, 
  const std::string& testName) const 
{
  bool rtn = true;
  if (err > spec.errorTol())
  {
    Out::root() << testName << " test FAILED: err=" << err << ", tol = " 
                << spec.errorTol() << std::endl;
    rtn = false;
  }
  else if (err > spec.warningTol())
  {
    Out::root() << "WARNING: " << testName << " test err="
                << err << " could not beat tol = " 
                << spec.warningTol() << std::endl;
  }
  else
  {
    Out::root() << "test " << testName << " PASSED with tol=" << spec.errorTol() << std::endl;
  }
  return rtn;
}
  
}
#endif
