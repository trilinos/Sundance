/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DISCRETESPACE_H
#define SUNDANCE_DISCRETESPACE_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceDOFMapBase.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFObjectWithVerbosity.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace TSFExtended;
  using namespace SundanceStdMesh;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Internal;

  /** 
   *
   */
  class DiscreteSpace : public ObjectWithVerbosity<DiscreteSpace>
  {
  public:
    /** */
    DiscreteSpace(){;}
    /** */
    DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
                  const VectorType<double>& vecType);
    /** */
    DiscreteSpace(const Mesh& mesh, const Array<BasisFamily>& basis,
                  const VectorType<double>& vecType);

    /** */
    DiscreteSpace(const Mesh& mesh, const Array<BasisFamily>& basis,
                  const RefCountPtr<DOFMapBase>& map,
                  const VectorType<double>& vecType);

    /** */
    const RefCountPtr<DOFMapBase>& map() const {return map_;}

    /** return the number of functions */
    int nFunc() const {return basis_.size();}

    /** */
    const Array<BasisFamily>& basis() const {return basis_;}

    /** */
    Vector<double> createVector() const {return vecSpace_.createMember();}

    /** */
    VectorSpace<double> vecSpace() const {return vecSpace_;}

    /** */
    const Mesh& mesh() const {return mesh_;}
  private:
    /** */
    RefCountPtr<DOFMapBase> map_;

    /** */
    Mesh mesh_;

    /** */
    Array<BasisFamily> basis_;

    /** */
    VectorSpace<double> vecSpace_;

    /** */
    VectorType<double> vecType_;
  };

}



#endif
