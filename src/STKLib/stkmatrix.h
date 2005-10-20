#ifndef __STKMATRIX_H
#define __STKMATRIX_H

#include <cstddef>
//#include <cstdlib>
#include <stdlib.h>
#include <stdexcept>
#include <cblas.h>

#define CHECKSIZE


//	class matrix_error : public std::logic_error {};/
//	class matrix_sizes_error : public matrix_error {};

/// defines a storage type
typedef enum
{
	STORAGE_UNDEFINED = 0,
	STORAGE_REGULAR,
	STORAGE_TRANSPOSED
} storage_type;


	/**
	 *  @brief Provides a matrix abstraction class
	 *
	 *  This class provides a way to work with matrices in STK.
	 *  It encapsulates basic operations and memory optimizations.
	 *
	 */
	template<typename _ElemT>
		class basic_stkmatrix
		{
		public:
			/// defines a type of this
			typedef basic_stkmatrix<_ElemT>    __this_type;


		protected:
			/// keeps info about data layout in the memory
			storage_type storageType;


			//@{
			/// these atributes store the real matrix size as it is stored in memory
			/// including memalignment
			size_t  M_rows;       ///< Number of rows
			size_t  M_cols;       ///< Number of columns
			size_t  M_realCols;   ///< true number of columns for the internal matrix.
			                      ///< This number may differ from M_cols as memory
														///< alignment might be used
			size_t  M_size;       ///< Total size of data block in bytes
			size_t  M_skip;       ///< Bytes to skip (memalign...)

			size_t  T_rows;       ///< Real number of rows (available to the user)
			size_t  T_cols;       ///< Real number of columns
			//@}

			/// data memory area
			_ElemT *  data;


		public:
			// Constructors

			/// Empty constructor
			basic_stkmatrix<_ElemT> (): storageType(STORAGE_UNDEFINED) {}

			/// Copy constructor
			basic_stkmatrix<_ElemT> (const __this_type & t);

			/// Basic constructor
			basic_stkmatrix<_ElemT> (const size_t r,
															 const size_t c,
															 const storage_type st = STORAGE_REGULAR);


			/// Destructor
			~basic_stkmatrix<_ElemT> ();


			/// Returns number of rows in the matrix
			virtual size_t
			rows() const
			{
				return T_rows;
			}

			/// Returns number of columns in the matrix
			virtual size_t
			cols() const
			{
				return T_cols;
			}


			/// Returns the way the matrix is stored in memory
			storage_type
			storage() const
			{
				return storageType;
			}


			/**
			 *  @brief Adds values of matrix a to this element by element
			 *  @param a Matrix to be added to this
			 *  @return Reference to this
			 */
			__this_type &
			operator += (const __this_type & a);

			/**
			 *  @brief Performs vector multiplication on a and b and and adds the
			 *         result to this (elem by elem)
			 */
			__this_type &
			AddVectorMult(const __this_type & a, const __this_type & b);


			/**
			 *  @brief Gives access to the matrix memory area
			 *  @return pointer to the first field
			 */
			virtual _ElemT *
			operator () () const {return data;};

			/**
			 *  @brief Gives access to a specified matrix row
			 *  @return pointer to the first field of the row
			 */
			virtual _ElemT *
			row (const size_t r)
			{
				return data + (r * M_realCols);
			}

			/**
			 *  @brief Gives access to matrix elements (row, col)
			 *  @return pointer to the desired field
			 */
			virtual _ElemT *
			operator () (const size_t r, const size_t c);

			void printOut();

		}; // class basic_stkmatrix



	/**
	 *  @brief Provides a window matrix abstraction class
	 *
	 *  This class provides a way to work with matrix cutouts in STK.
	 *  It encapsulates basic operations and memory optimizations.
	 *
	 */
	template<typename _ElemT>
		class basic_stkwinmatrix : public basic_stkmatrix<_ElemT>
		{
		private:
			/// points to the original begining of the data array
			/// The data atribute points now to the begining of the window
			_ElemT * subdata;

			size_t T_rows;    ///< Real number of window rows
			size_t T_cols;    ///< Real number of window columns
			size_t T_rowOff;  ///< First row of the window
			size_t T_colOff;  ///< First column of the window

			size_t M_skip;    ///< Columns to skip in the memory



		public:
			/// defines a type of this
			typedef basic_stkwinmatrix<_ElemT>    __this_type;


			/// Empty constructor
			basic_stkwinmatrix() : basic_stkmatrix<_ElemT>() {};

			/// Copy constructor
			basic_stkwinmatrix<_ElemT> (const __this_type & t);

			/// Basic constructor
			basic_stkwinmatrix<_ElemT> (const size_t r,
																	const size_t c,
																	const storage_type st = STORAGE_REGULAR):
				basic_stkmatrix<_ElemT>(r, c, st),
				T_rows(r),   T_cols(c),
				T_rowOff(0), T_colOff(0), M_skip(basic_stkmatrix<_ElemT>::M_skip)
			{
				subdata = basic_stkmatrix<_ElemT>::data;
			}


			void
			setSize(const size_t _ro,
							const size_t _r,
							const size_t _co=0,
							const size_t _c=basic_stkmatrix<_ElemT>::T_cols);

			/**
			 *  @brief Resets the window to the default (full) size
			 */
			void
			reset ()
			{
				subdata  = basic_stkmatrix<_ElemT>::data;

				T_rows   = basic_stkmatrix<_ElemT>::T_rows;
				T_cols   = basic_stkmatrix<_ElemT>::T_cols;
				T_rowOff = 0;
				T_colOff = 0;
				M_skip   = basic_stkmatrix<_ElemT>::M_skip;
			}

			/**
			 *  @brief Gives access to the matrix memory area
			 *  @return pointer to the first field
			 */
			virtual _ElemT *
			operator () () const {return subdata;};


			/**
			 *  @brief Gives access to a specified matrix row
			 *  @return pointer to the first field of the row
			 */
			//virtual _ElemT *
			//row (const size_t r);

			/// this operator will print out the whole matrix
			//friend ostream& operator << (ostream & out, __this_type & m);
			//void printOut();




		};

// we need to include the implementation
#include "stkmatrix.tcc"



#endif //#ifndef __STKMATRIX_H
