/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DISCRETEFUNCTION_H
#define SUNDANCE_DISCRETEFUNCTION_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceFuncWithBasis.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "TSFVector.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Internal;


  /** 
   *
   */
  class DiscreteFunction : public DiscreteFunctionStub,
                           public FuncWithBasis
  {
  public:
    /** */
    DiscreteFunction(const DiscreteSpace& space, const string& name="");

    /** */
    DiscreteFunction(const DiscreteSpace& space, const Vector<double>& vec, 
                     const string& name="");

    /** */
    DiscreteFunction(const DiscreteSpace& space, const double& constantValue,
                     const string& name="");
   

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** virtual destructor */
    virtual ~DiscreteFunction() {;}

    /* boilerplate */
    GET_RCP(ExprBase);

    /** */
    const Vector<double>& vector() const {return vector_;}

    /** */
    void setVector(const Vector<double>& vec) {vector_ = vec;}

    /** */
    const DiscreteSpace& discreteSpace() const {return space_;}

    /** */
    const Mesh& mesh() const {return space_.mesh();}

    /** */
    const RefCountPtr<DOFMapBase>& map() const {return space_.map();}

    /** */
    void getLocalValues(int cellDim, 
                        const Array<int>& cellLID,
                        Array<double>& localValues) const ;


  private:

    DiscreteSpace space_;

    Vector<double> vector_;

#endif /* DOXYGEN_DEVELOPER_ONLY */
  };

}



#endif
