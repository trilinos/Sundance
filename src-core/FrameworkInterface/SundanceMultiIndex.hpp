/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MULTIINDEX_H
#define SUNDANCE_MULTIINDEX_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_XMLObject.hpp"
#include <string>
#include <stdexcept>

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using std::runtime_error;
  using namespace Teuchos;
  using std::string;
  using std::ostream;

  namespace Internal
    {
      /**
       * An integer vector representing a multivariate derivative.
       */

      class MultiIndex
        {
        public:
          /** constructs D(0,0,0) */
          MultiIndex();
          /** constructs a multiindex D(x,y,z) */
          MultiIndex(int x, int y, int z);

          /** */
          bool operator==(const MultiIndex& other) const ;

          /** */
          bool operator<(const MultiIndex& other) const ;

          /** */
          const int& operator[](int i) const {return m_[i];}

          /** */
          int& operator[](int i) {return m_[i];}

          /** */
          MultiIndex operator+(const MultiIndex& other) const ;

          /** Return a list of all multiindices of lower order than this one. */
          Array<MultiIndex> getLowerMultiIndices() const ;

          /** */
          virtual string toString() const ;

          /** */
          XMLObject toXML() const ;

          /** */
          virtual int hashCode() const ;

          /** */
          int order() const ;

          /** */
          int firstOrderDirection() const ;

          /** */
          static int maxDim() {return 3;}

          /** */
          string coordForm() const ;
        private:
          Array<int> m_;
        };

    }
}

namespace Teuchos
{
  using std::string;
  /** \relates SundanceCore::Internal::MultiIndex */
  inline int hashCode(const SundanceCore::Internal::MultiIndex& h)
    {return h.hashCode();}

  /** \relates SundanceCore::Internal::MultiIndex */
  inline string toString(const SundanceCore::Internal::MultiIndex& h)
    {return h.toString();}
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
