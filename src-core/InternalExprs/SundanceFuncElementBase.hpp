/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FUNCELEMENTBASE_H
#define SUNDANCE_FUNCELEMENTBASE_H


#include "SundanceDefs.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceSet.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  namespace Internal
    {
      /** 
       * FuncElementBase defines the interface for scalar-valued elements
       * of Sundance functions. At the user level, Sundance functions can be
       * list (e.g, vector or tensor) valued; internally, however, compound
       * expressions use only scalar functions deriving from the 
       * FuncElementBase class. 
       *
       * <h4> Function ID </h4>
       *
       * Every function element has a unique integer 
       * identifier, or <b> funcID, </b>
       * which is used internally as a shorthand for the function. The funcID
       * is assigned at construction using the private static method
       * nextID(), and because it is private with no mutators, can never
       * be changed. The nextID() method is responsible for guaranteeing that 
       * funcIDs are unique.
       */
      class FuncElementBase : virtual public ScalarExpr
        {
        public:
          /** */
          FuncElementBase(const string& rootName,
                          const string& suffix);
          /** */
          FuncElementBase(const string& rootName);

          /** virtual destructor */
          virtual ~FuncElementBase() {;}

          /** Return an integer ID which uniquely identifies this function */
          int funcID() const {return id_;}

          /** Append to the set of func IDs present in 
           * this expression. */
          virtual void accumulateFuncSet(Set<int>& funcIDs,
                                         const Set<int>& activeFuncs) const ;

          /** Return the name of this function */
          const string& name() const {return name_;}

          /** Return the root name of this function */
          const string& rootName() const {return rootName_;}

          /** Return the root name of this function */
          const string& suffix() const {return suffix_;}

          /** Write self in text form */
          virtual ostream& toText(ostream& os, bool paren) const ;

          /** Write self in Latex form */
          virtual ostream& toLatex(ostream& os, bool paren) const ;

        protected:
          /** Determine whether this function is in the given active set */
          bool isInActiveSet(const Set<MultiSet<int> >& activeFuncIDs) const ;
          
        private:

          string name_;

          string rootName_;

          string suffix_;

          int id_;

          static int& nextID() {static int rtn = 0; return rtn;}

        };
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
