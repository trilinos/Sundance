/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_REGIONQUADCOMBO_H
#define SUNDANCE_REGIONQUADCOMBO_H


#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_Utils.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceOrderedHandle.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY


using SundanceUtils::Map;
namespace SundanceCore
{
  using namespace Teuchos;
  using namespace SundanceUtils;
  using std::string;


  namespace Internal
    {

      /** */
      typedef OrderedPair<OrderedHandle<CellFilterStub>,
                          OrderedHandle<QuadratureFamilyStub> > RegPair;
      /** 
       * Expressions may appear in more than one subregions of a problem,
       * for instance in an internal domain and also on a boundary. On
       * those different subregions, a given expression might be subject
       * to different sets of functional derivatives; thus, different
       * evaluation regions might have different sparsity patterns. 
       * It is therefore necessary to build and store
       * sparsity information on a region-by-region basis. 
       *
       * Class RegionQuadCombo is used as an identifier for regions. The
       * only thing it needs to do is to be useable as a key in a STL map.
       */
      class RegionQuadCombo 
        {
        public:
          /** */
          RegionQuadCombo();
          /** */
          RegionQuadCombo(const RefCountPtr<CellFilterStub>& domain,
                          const RefCountPtr<QuadratureFamilyStub>& quad);

          /** */
          inline bool operator==(const RegionQuadCombo& other) const
            {return id_==other.id_;}

          /** */
          string toString() const ;

          /** */
          bool operator<(const RegionQuadCombo& other) const
            {return id_ < other.id_;}

          /** */
          const RefCountPtr<CellFilterStub>& domain() const {return domain_;}

          /** */
          const RefCountPtr<QuadratureFamilyStub>& quad() const 
          {return quad_;}

        private:

          /** */
          int id_;

          /** */
          RefCountPtr<CellFilterStub> domain_;

          /** */
          RefCountPtr<QuadratureFamilyStub> quad_;
          
          /** */
          static int getID(const RefCountPtr<CellFilterStub>& domain,
                           const RefCountPtr<QuadratureFamilyStub>& quad);

          /** */
          static int topID() {static int rtn=0; return rtn++;}

          /** */
          static Map<RegPair, int>& domainAndQuadToIDMap() ;
        };

    }
}

namespace std
{
  /** \relates SundanceCore::Internal::RegionQuadCombo*/
  inline ostream& operator<<(ostream& os, 
                             const SundanceCore::Internal::RegionQuadCombo& c)
  {
    os << c.toString();
    return os;
  }
}

namespace Teuchos
{
  using std::string;

  /** \relates SundanceCore::Internal::RegionQuadCombo */
  inline string toString(const SundanceCore::Internal::RegionQuadCombo& h)
    {return h.toString();}

}





#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
