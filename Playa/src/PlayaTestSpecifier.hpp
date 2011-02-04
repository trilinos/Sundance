/* @HEADER@ */
//   
 /* @HEADER@ */


#ifndef PLAYA_TESTSPECIFIER_HPP
#define PLAYA_TESTSPECIFIER_HPP

namespace Playa
{
  /** */
  template <class Scalar>
  class TestSpecifier
  {
  public:
    /** \brief Local typedef for promoted scalar magnitude */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** 
     * 
     */
    TestSpecifier(bool doTest, 
                  const ScalarMag& warningTol,
                  const ScalarMag& errorTol)
      : doTest_(doTest),
        warningTol_(warningTol),
        errorTol_(errorTol)
    {;}

    /** */
    bool doTest() const {return doTest_;}

    /** */
    const ScalarMag& warningTol() const {return warningTol_;}

    /** */
    const ScalarMag& errorTol() const {return errorTol_;}

  private:

    bool doTest_;

    double warningTol_;

    double errorTol_;
  };
}
#endif
