/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_QUADRATUREFAMILYBASE_H
#define SUNDANCE_QUADRATUREFAMILYBASE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceCellType.hpp"
#include "SundancePoint.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace Teuchos;

  namespace Internal
  {
    /** 
     * QuadratureFamilyBase extends QuadratureFamilyStub to provide
     * an interface for getting quadrature points for a given cell type. 
     */
    class QuadratureFamilyBase : public QuadratureFamilyStub 
    {
    public:
      /** */
      QuadratureFamilyBase(int order) : QuadratureFamilyStub(order) {;}

      /** */
      virtual ~QuadratureFamilyBase(){;}

      /** Get the quadrature points and weights for the given cell type */
      void getPoints(const CellType& cellType, 
                     Array<Point>& quadPoints,
                     Array<double>& quadWeights) const ;
    protected:

      /** compute a rule for the reference line cell */
      virtual void getLineRule(Array<Point>& quadPoints,
                               Array<double>& quadWeights) const ;

      /** compute a rule for the reference triangle cell */
      virtual void getTriangleRule(Array<Point>& quadPoints,
                                   Array<double>& quadWeights) const ;

      /** compute a rule for the reference quad cell */
      virtual void getQuadRule(Array<Point>& quadPoints,
                               Array<double>& quadWeights) const ;

      /** compute a rule for the reference tet cell */
      virtual void getTetRule(Array<Point>& quadPoints,
                              Array<double>& quadWeights) const ;


      /** compute a rule for the reference brick cell */
      virtual void getBrickRule(Array<Point>& quadPoints,
                                Array<double>& quadWeights) const ;
    private:
    };
  }
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif

