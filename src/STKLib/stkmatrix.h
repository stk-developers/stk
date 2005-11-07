#ifndef __STKMATRIX_H
#define __STKMATRIX_H

#include <cstddef>
//#include <cstdlib>
#include <stdlib.h>
#include <stdexcept>
#include <cblas.h>

#include <iostream>

#define CHECKSIZE


namespace STK
{

  //  class matrix_error : public std::logic_error {};/
  //  class matrix_sizes_error : public matrix_error {};

  /// defines a storage type
  typedef enum
  {
    STORAGE_UNDEFINED = 0,
    STORAGE_REGULAR,
    STORAGE_TRANSPOSED
  } storage_type;


    // declare the class so the header knows about it
    template<typename _ElemT>
      class Matrix;


    // we need to declare the friend << operator here
    template<typename _ElemT>
      std::ostream & operator << (std::ostream & out, const Matrix<_ElemT> & m);


    /**
     *  @brief Provides a matrix abstraction class
     *
     *  This class provides a way to work with matrices in STK.
     *  It encapsulates basic operations and memory optimizations.
     *
     */
    template<typename _ElemT>
      class Matrix
      {
      public:
        /// defines a type of this
        typedef Matrix<_ElemT>    __this_type;


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
        Matrix<_ElemT> (): storageType(STORAGE_UNDEFINED) {}

        /// Copy constructor
        Matrix<_ElemT> (const __this_type & t);

        /// Basic constructor
        Matrix<_ElemT> (const size_t r,
                                 const size_t c,
                                 const storage_type st = STORAGE_REGULAR);

        /// Destructor
        ~Matrix<_ElemT> ();


        /// Returns number of rows in the matrix
        size_t
        rows() const
        {
          return T_rows;
        }

        /// Returns number of columns in the matrix
        size_t
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
        _ElemT *
        operator () () {return data;};

        /**
         *  @brief Gives access to a specified matrix row
         *  @return pointer to the first field of the row
         */
        _ElemT *
        row (const size_t r)
        {
          return data + (r * M_realCols);
        }

        /**
         *  @brief Gives access to matrix elements (row, col)
         *  @return pointer to the desired field
         */
        _ElemT *
        operator () (const size_t r, const size_t c);


        friend std::ostream & operator << <> (std::ostream & out, const __this_type & m);

        void printOut();

      }; // class Matrix



    /**
     *  @brief Provides a window matrix abstraction class
     *
     *  This class provides a way to work with matrix cutouts in STK.
     *  It encapsulates basic operations and memory optimizations.
     *
     */
    template<typename _ElemT>
      class WindowMatrix : public Matrix<_ElemT>
      {
      protected:
        /// points to the original begining of the data array
        /// The data atribute points now to the begining of the window
        _ElemT * orig_data;

        //@{
        /// these atributes store the real matrix size as it is stored in memory
        /// including memalignment
        size_t  orig_M_rows;       ///< Number of rows
        size_t  orig_M_cols;       ///< Number of columns
        size_t  orig_M_realCols;   ///< true number of columns for the internal matrix.
                                   ///< This number may differ from M_cols as memory
                                   ///< alignment might be used
        size_t  orig_M_size;       ///< Total size of data block in bytes
        size_t  orig_M_skip;       ///< Bytes to skip (memalign...)

        size_t  orig_T_rows;       ///< Original number of window rows
        size_t  orig_T_cols;       ///< Original number of window columns

        size_t T_rowOff;           ///< First row of the window
        size_t T_colOff;           ///< First column of the window
        //@}


      public:
        /// defines a type of this
        typedef WindowMatrix<_ElemT>    __this_type;



        /// Empty constructor
        WindowMatrix() : Matrix<_ElemT>() {};

        /// Copy constructor
        WindowMatrix<_ElemT> (const __this_type & t);

        /// Basic constructor
        WindowMatrix<_ElemT> (const size_t r,
                                    const size_t c,
                                    const storage_type st = STORAGE_REGULAR):
          Matrix<_ElemT>(r, c, st), // create the base class
          orig_M_rows    (Matrix<_ElemT>::M_rows),
          orig_M_cols    (Matrix<_ElemT>::M_cols),
          orig_M_realCols(Matrix<_ElemT>::M_realCols),
          orig_M_size    (Matrix<_ElemT>::M_size),
          orig_M_skip    (Matrix<_ElemT>::M_skip),
          orig_T_rows    (Matrix<_ElemT>::T_rows),                    // copy the original values
          orig_T_cols    (Matrix<_ElemT>::T_cols),
          T_rowOff(0), T_colOff(0)          // set the offset

        {
          orig_data = Matrix<_ElemT>::data;
        }

        /// The destructor
        ~WindowMatrix<_ElemT>()
        {
          Matrix<_ElemT>::data = this->orig_data;
        }

        /// sets the size of the window (whole rows)
        void
        setSize(const size_t _ro,
                const size_t _r);


        /// sets the size of the window
        void
        setSize(const size_t _ro,
                const size_t _r,
                const size_t _co,
                const size_t _c);

        /**
         *  @brief Resets the window to the default (full) size
         */
        void
        reset ()
        {
          Matrix<_ElemT>::data      =  orig_data;
          Matrix<_ElemT>::T_rows    =  orig_T_rows;                        // copy the original values
          Matrix<_ElemT>::T_cols    =  orig_T_cols;
          Matrix<_ElemT>::M_rows    =  orig_M_rows;
          Matrix<_ElemT>::M_cols    =  orig_M_cols;
          Matrix<_ElemT>::M_realCols=  orig_M_realCols;
          Matrix<_ElemT>::M_size    =  orig_M_size;
          Matrix<_ElemT>::M_skip    =  orig_M_skip;

          T_rowOff = 0;
          T_colOff = 0;
        }
      };

} // namespace STK
// we need to include the implementation
#include "stkmatrix.tcc"

#endif //#ifndef __STKMATRIX_H
