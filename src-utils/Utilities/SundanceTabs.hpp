/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TABS_H
#define SUNDANCE_TABS_H

#include "SundanceDefs.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceUtils
{
  /**
   * Tabbing utility for output. Constructing a new Tabs object automatically
   * increments the number of tabs to be written. When the Tabs object goes out
   * of scope, the original tabs level is restored. 
   *
   * The tab size and character can be specified through the setTabSize() and
   * setTabChar() methods, for example,
   * \code
   * Tabs::setTabChar('*');
   * Tabs::setTabSize(4);
   * \endcode
   * The tab character can be set on an object-by-object basis
   * through a constructor argument.
   *
   * By default, a header giving the depth of tabs is written to each line; this
   * can simplify scanning by eye for when a given tab level is reached. 
   * This header can be turned off by calling
   * \code
   * Tabs::showDepth() = false;
   * \endcode
   * 
   * Example: the code
   * \code
   * void f()
   * {
   *   Tabs tab;
   *   cout << tab << "in f()" << endl;
   *   g();
   *   cout << tab << "leaving f()" << endl;
   * }
   *
   * void g()
   * {
   *   Tabs tab0;
   *   cout << tab0 << "in g()" << endl;
   *   for (int i=0; i<3; i++)
   *     {
   *       Tabs tab1('-');
   *       cout << tab1 << "i=" << i << endl;
   *     }
   *   cout << tab0 << "leaving g()" << endl;
   * }
   * \endcode
   * writes the following output 
   * \code
   * [0]  in f()
   * [1]    in g()
   * [2]------i=0
   * [2]------i=1
   * [2]------i=2
   * [1]    leaving g()
   * [0]  leaving f()
   * \endcode 
   */
  class Tabs
    {
    public:
      /** Constructor increments tab level */
      Tabs(char c = tabChar());

      /** Destructor decrements tab level */
      ~Tabs();

      /** 
       * Print to stream. This method is usually not called directly, as
       * tabs will usually be written with the insertion operator
       */
      void print(ostream& os) const ;

      /** Change the tab size. Default is 2.  */
      static void setTabSize(int ts) {tabSize() = ts;}

      /** Change the default tab character. Default is an empty space. */
      static void setTabChar(char tc) {tabChar() = tc;}

      /** Indicate whether to print the tab depth as a header for each line. */
      static bool& showDepth() {static bool rtn = true; return rtn;}

    private:
      /** */
      static int& tabLevel() {static int rtn = 0; return rtn;}

      /** */
      static int& tabSize() {static int rtn = 2; return rtn;}

      /** */
      static char& tabChar() {static char rtn = ' '; return rtn;}

      char c_;
    };
}

namespace std
{
  /** \relates Tabs stream insertion operator for tab */
  inline ostream& operator<<(ostream& os, const SundanceUtils::Tabs& t)
    {
      t.print(os);
      return os;
    }
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif
