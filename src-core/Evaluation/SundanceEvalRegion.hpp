/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALREGION_H
#define SUNDANCE_EVALREGION_H


#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_Utils.hpp"
#include "SundanceCellFilterBase.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceOrderedHandle.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY



namespace SundanceCore
{
  using namespace Teuchos;
  using namespace SundanceUtils;
  using std::string;


  namespace FrameworkInterface
    {

      /** */
      typedef OrderedPair<OrderedHandle<CellFilterBase>,
                          OrderedHandle<QuadratureFamilyBase> > RegPair;
      /** 
       * Expressions may appear in more than one subregions of a problem,
       * for instance in an internal domain and also on a boundary. On
       * those different subregions, a given expression might be subject
       * to different sets of functional derivatives; thus, different
       * evaluation regions might have different sparsity patterns. 
       * It is therefore necessary to build and store
       * sparsity information on a region-by-region basis. 
       *
       * Class EvalRegion is used as an identifier for regions. The
       * only thing it needs to do is to be useable as a key in a STL map.
       */
      class EvalRegion
        {
        public:
          /** */
          EvalRegion();
          /** */
          EvalRegion(const RefCountPtr<CellFilterBase>& domain,
                     const RefCountPtr<QuadratureFamilyBase>& quad);

          /** */
          inline bool operator==(const EvalRegion& other) const
            {return id_==other.id_;}

          /** */
          string toString() const ;

          /** */
          bool operator<(const EvalRegion& other) const
            {return id_ < other.id_;}

        private:

          /** */
          int id_;

          /** */
          RefCountPtr<CellFilterBase> domain_;

          /** */
          RefCountPtr<QuadratureFamilyBase> quad_;
          
          /** */
          static int getID(const RefCountPtr<CellFilterBase>& domain,
                           const RefCountPtr<QuadratureFamilyBase>& quad);

          /** */
          static int topID() {static int rtn=0; return rtn++;}

          /** */
          static Map<RegPair, int>& domainAndQuadToIDMap() ;
        };

    }
}

namespace std
{
  /** \relates SundanceCore::FrameworkInterface::EvalRegion*/
  inline ostream& operator<<(ostream& os, 
                             const SundanceCore::FrameworkInterface::EvalRegion& c)
  {
    os << c.toString();
    return os;
  }
}

namespace Teuchos
{
  using std::string;

  /** \relates SundanceCore::FrameworkInterface::EvalRegion */
  inline string toString(const SundanceCore::FrameworkInterface::EvalRegion& h)
    {return h.toString();}

}





#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
