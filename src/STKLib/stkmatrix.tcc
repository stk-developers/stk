#include <cstdlib>
#include <cblas.h>


//******************************************************************************
	// Basic constructor
	template<typename _ElemT>
	basic_stkmatrix<_ElemT>::
	basic_stkmatrix<_ElemT> (const size_t       _r,
													 const size_t       _c,
													 const storage_type _st)
	{
		// initialize some helping vars
		size_t  _rows;
		size_t  _cols;
		size_t  _skip;
		size_t  _realcols;
		size_t  _size;
		void ** _data;

		// at this moment, nothing is clear
		storageType = STORAGE_UNDEFINED;

		// if we store transposed, we swap the rows and cols
		// we assume that the user does not specify STORAGE_UNDEFINED
		if (_st == STORAGE_TRANSPOSED)
		{ _rows = _c; _cols = _r; }
		else
		{ _rows = _r; _cols = _c; }


		// compute the size of skip and real cols
		_skip     = ((16 / sizeof(_ElemT)) - _cols % (16 / sizeof(_ElemT))) % (16 / sizeof(_ElemT));
		_realcols = _cols + _skip;
		_size     = _rows * _realcols * sizeof(_ElemT);


		// allocate the memory and set the right dimensions and parameters
		if (!posix_memalign(_data, 16, _size))
		{
			this->data  = static_cast<_ElemT *> (*_data);

			storageType = _st;
			M_rows      = _rows;
			M_cols      = _cols;
			M_realCols  = _realcols;
			M_skip      = _skip;
			M_size      = _size;

			T_rows      = _r;
			T_cols      = _c;
		}
		else
		{
			throw std::bad_alloc();
		}
	}


//******************************************************************************
	// Copy constructor
	template<typename _ElemT>
	basic_stkmatrix<_ElemT>::
	basic_stkmatrix<_ElemT> (const __this_type & t)
	{
		// we need to be sure that the storage type is defined
		if (t.st != STORAGE_UNDEFINED)
		{
			void ** _data;

			// first allocate the memory
			// this->data = new char[t.M_size];
			// allocate the memory and set the right dimensions and parameters
			if (!posix_memalign(_data, 16, t.M_size))
			{
				this->data  = static_cast<_ElemT *> (*_data);

				// copy the memory block
				memcpy(this->data, t.data, t.M_size);

				// set the parameters
				this->storageType = t.storageType;
				this->M_rows      = t.M_rows;
				this->M_cols      = t.M_cols;
				this->M_realCols  = t.M_realCols;
				this->M_skip      = t.M_skip;
				this->M_size      = t.M_size;
			}
			else
			{
				// throw bad_alloc exception if failure
				throw std::bad_alloc();
			}
		}
	}


//******************************************************************************
	template<typename _ElemT>
	basic_stkmatrix<_ElemT> &
	basic_stkmatrix<_ElemT>::
	operator += (const basic_stkmatrix<_ElemT> & a)
	{
#ifdef CHECKSIZE
		// the sizes must match, else throw an exception
		if (a.M_size != this->M_size)
			throw std::logic_error("Sizes of matrices must be identical");
#endif

		// return a reference to this
		return *this;
	}


//******************************************************************************
	// The destructor
	template<typename _ElemT>
	basic_stkmatrix<_ElemT>::
	~basic_stkmatrix<_ElemT> ()
	{
		// we need to free the data block if it was defined
		if (storageType != STORAGE_UNDEFINED)
		{
			free(data);
		}
	}

//******************************************************************************
	template<typename _ElemT>
	_ElemT *
	basic_stkmatrix<_ElemT>::
	operator ()(const size_t r, const size_t c)
	{
		//TODO
		//
	}



//******************************************************************************
	// Copy constructor
	template<typename _ElemT>
		basic_stkwinmatrix<_ElemT>::
		basic_stkwinmatrix<_ElemT> (const __this_type & t) :
			basic_stkmatrix<_ElemT> (t)
		{
			// we need to be sure that the storage type is defined
			if (t.st != STORAGE_UNDEFINED)
			{
				// the data should be copied by the parent constructor
				void ** _data;

				// first allocate the memory
				// this->data = new char[t.M_size];
				// allocate the memory and set the right dimensions and parameters
				if (!posix_memalign(_data, 16, t.M_size))
				{
					this->data  = static_cast<_ElemT *> (*_data);

					// copy the memory block
					memcpy(this->data, t.data, t.M_size);

					// set the parameters
					this->storageType = t.storageType;
					this->M_rows      = t.M_rows;
					this->M_cols      = t.M_cols;
					this->M_realCols  = t.M_realCols;
					this->M_skip      = t.M_skip;
					this->M_size      = t.M_size;
				}
				else
				{
					// throw bad_alloc exception if failure
					throw std::bad_alloc();
				}
			}
		}


//******************************************************************************
	// Copy constructor
	template<typename _ElemT>
		void
		basic_stkwinmatrix<_ElemT>::
		setSize(const size_t _ro,
						const size_t _r,
						const size_t _co,
						const size_t _c)
		{
			if (basic_stkmatrix<_ElemT>::storageType == STORAGE_REGULAR)
			{
				// point to the begining of window
				subdata = basic_stkmatrix<_ElemT>::data +
									basic_stkmatrix<_ElemT>::M_realCols * _ro + _co;

				this->T_rows = _r;
				this->T_cols = _c;
			}
			else if (basic_stkmatrix<_ElemT>::storageType == STORAGE_TRANSPOSED)
			{
				// point to the begining of window
				subdata = basic_stkmatrix<_ElemT>::data +
									basic_stkmatrix<_ElemT>::M_realCols * _co + _ro;

				this->T_rows = _r;
				this->T_cols = _c;
			}
		}

using namespace std;

	template<typename _ElemT>
		void
		basic_stkmatrix<_ElemT>::
		printOut()
		//ostream& operator << (ostream& out, basic_stkmatrix<_ElemT> & m)
		{
			std::cerr << "Starting" << endl;
			size_t    rcount;
			size_t    ccount;

			_ElemT  * data    = this->data;

			// go through all rows
			for (rcount = 0; rcount < this->M_rows; rcount++)
			{
				// go through all columns
				for (ccount = 0; ccount < this->M_cols; ccount++)
				{
					std::cout << static_cast<_ElemT> (*data) << " ";
					data ++;
				}

				std::cout << std::endl;
			}
		}

