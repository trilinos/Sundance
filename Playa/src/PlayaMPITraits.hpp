// @HEADER
// @HEADER

#ifndef PLAYA_MPITRAITS_H
#define PLAYA_MPITRAITS_H

/*! \file PlayaMPITraits.hpp
 * \brief Declaration of a templated traits class for binding MPI types to
 * C++ types. This is for use with the MPIComm class and is supposed to compile
 * rgeardless of whether MPI has been enabled. If you need to convert directly to 
 * MPI types (e.g., MPI_INT), please refer to Teuchos_MPIRawTraits.hpp.
*/

#include "PlayaMPIComm.hpp"

namespace Playa
{
	using std::string;

	/** \ingroup MPI 
	 * \brief Templated traits class that binds MPI types to C++ types
	 * \note Template specializations exist for datatypes: <tt>char</tt>,
		<tt>int</tt>, <tt>float</tt>, and <tt>double</tt>.
	 */
	template <class T> class MPITraits
		{
		public:
			/** \brief Return the MPI data type of the template argument */
			static MPIDataType type();
		};

	/** \ingroup MPI 
	 */
	template <> class MPITraits<int>
		{
		public:
			/** return the MPI data type of the template argument */
			static MPIDataType type() {return MPIDataType::intType();}
		};
	
	/** \ingroup MPI 
	 */
	template <> class MPITraits<float>
		{
		public:
			/** return the MPI data type of the template argument */
			static MPIDataType type() {return MPIDataType::floatType();}
		};
	
	/** \ingroup MPI 
	 */
	template <> class MPITraits<double>
		{
		public:
			/** return the MPI data type of the template argument */
			static MPIDataType type() {return MPIDataType::doubleType();}
		};
	
	/** \ingroup MPI 
	 */
	template <> class MPITraits<char>
		{
		public:
			/** return the MPI data type of the template argument */
			static MPIDataType type() {return MPIDataType::charType();}
		};

	
} // namespace Playa

#endif
