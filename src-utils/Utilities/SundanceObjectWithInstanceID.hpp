/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_OBJECTWITHINSTANCEID_H
#define SUNDANCE_OBJECTWITHINSTANCEID_H

#include "SundanceDefs.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceUtils
{
    /**
     * ObjectWithInstanceID provides a common method for the
     * generation of instance-specific ID numbers. Subclasses will inherit
     * the id() method, and instances of those subclasses will be
     * given a unique ID at construction time.
     * 
     * <h4> Design note: </h4> By templating on the derived type, 
     * we can give each derived type its own sequence of ID numbers.
     */
    template <class T>
    class ObjectWithInstanceID
    {
    public:
      /** Empty ctor will assign ID at construction time */
      ObjectWithInstanceID() : id_(nextID()) {;} 

      /** Return this object's ID number */
      int id() const {return id_;}

    private:
      /** Generate the next ID in the sequence */
      static int& nextID() {static int rtn=0; rtn++; return rtn;}

      /** */
      int id_;
    };

}


#endif  /* DOXYGEN_DEVELOPER_ONLY */   
#endif
