/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ARRAYOFTUPLES_H
#define SUNDANCE_ARRAYOFTUPLES_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "Teuchos_Array.hpp"

namespace SundanceUtils
{
  using namespace Teuchos;

  /** 
   * Class ArrayOfTuples packs an heterogeneous array of tuples into 
   * a single 1D array. 
   */
  template<class T>
    class ArrayOfTuples
    {
    public:
      /** Empty ctor */
      ArrayOfTuples();

      /** Constructor specifying the size of each tuple, but not the 
       * number of tuples */
      ArrayOfTuples(int tupleSize);

      /** Constructor specifying both the size and number of the tuples */
      ArrayOfTuples(int numTuples, int tupleSize);

      /** Returns the number of tuples */
      int length() const {return numTuples_;}

      /** Returns the size of the tuples */
      int tupleSize() const {return tupleSize_;}

      /** Change the number of tuples */
      void resize(int newSize) {data_.resize(newSize*tupleSize_);}

      /** Reserve memory for a number of tuples */
      void reserve(int newSize) {data_.reserve(newSize*tupleSize_);}

      /** Specify the size of the tuples */
      void setTupleSize(int tupleSize) {tupleSize_ = tupleSize;}

      /** Get the j-th entry in the i-th tuple */
      const T& value(int i, int j) const {return data_[i*tupleSize_+j];}

      /** Get the j-th entry in the i-th tuple */
      T& value(int i, int j) {return data_[i*tupleSize_+j];}

      /** Append a new tuple to the array */
      void append(const Array<T>& x);

      /** Append a new tuple to the array */
      void append(const T* x, int n);

    private:

      int numTuples_;
      int tupleSize_;

      Array<T> data_;

    };

  template<class T> inline ArrayOfTuples<T>::ArrayOfTuples()
    : numTuples_(0), tupleSize_(0), data_()
    {;}

  template<class T> inline ArrayOfTuples<T>::ArrayOfTuples(int tupleSize)
    : numTuples_(0), tupleSize_(tupleSize), data_()
    {;}

  template<class T> inline ArrayOfTuples<T>::ArrayOfTuples(int numTuples, int tupleSize)
    : numTuples_(numTuples), tupleSize_(tupleSize), data_(numTuples*tupleSize)
    {;}

  template<class T> inline void ArrayOfTuples<T>::append(const Array<T>& x)
    {
      for (int i=0; i<x.length(); i++)
        {
          data_.append(x[i]);
        }
      numTuples_++;
    }

  template<class T> inline void ArrayOfTuples<T>::append(const T* x, int n)
    {
      for (int i=0; i<n; i++)
        {
          data_.append(x[i]);
        }
      numTuples_++;
    }


}

#endif  /* DOXYGEN_DEVELOPER_ONLY */   

#endif
