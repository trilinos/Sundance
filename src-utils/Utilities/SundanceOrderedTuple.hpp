/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ORDEREDTUPLE_H
#define SUNDANCE_ORDEREDTUPLE_H

#include "SundanceDefs.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceUtils
{
  using namespace Teuchos;

  /** OrderedPair provides a means of lexigraphic comparison of a pair of
   * objects. The pair {a1, b1} is compared to {a2, b2} by first
   * comparing the most significant entries a1 and a2, and if they are
   * equal, comparing the least significant entries b1 and b2. */
  template<class A, class B>
    class OrderedPair
    {
    public:
      /** */
      OrderedPair(const A& a, const B& b)
        : a_(a), b_(b) {;}

      /** */
      inline bool operator<(const OrderedPair<A, B>& other) const
        {
          if ( a_ < other.a_ ) return true;
          if ( other.a_ < a_) return false;

          return b_ < other.b_;
        }

    private:
      const A a_;
      const B b_;
    };

  /** Lexigraphically-comparable triple of objects. */
  template<class A, class B, class C>
    class OrderedTriple : public OrderedPair<A, OrderedPair<B, C> >
    {
    public:
      /** */
      OrderedTriple(const A& a, const B& b, const C& c)
        : OrderedPair<A, OrderedPair<B, C> >(a, OrderedPair<B,C>(b,c))
        {;}
    };

  /** Lexigraphically-comparable quartet of objects. */
  template<class A, class B, class C, class D>
    class OrderedQuartet : public OrderedPair<A, OrderedTriple<B, C, D> >
    {
    public:
      /** */
      OrderedQuartet(const A& a, const B& b, const C& c, const D& d)
        : OrderedPair<A, OrderedTriple<B, C, D> >(a, OrderedTriple<B,C,D>(b,c,d))
        {;}
    };

}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
