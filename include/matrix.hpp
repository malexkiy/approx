#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

#include <vector>
#include <algorithm>

namespace nma
{
	template<typename T>
	class matrix
	{
	private:
		std::vector<T> _data;
		unsigned long _rows;
		unsigned long _cols;

	public:

		matrix(unsigned long rows, unsigned long cols);

		unsigned long rows() const;
		unsigned long cols() const;
		std::vector<T> data() const;

		static matrix<T> identity(unsigned long size);

		T& operator()(unsigned long row, unsigned long col) const;
		matrix<T> operator*(const matrix& other) const;
		matrix<T> transpose() const;
	};


	template<typename T>
	matrix<T>::matrix(unsigned long rows, unsigned long cols) :
		_rows(rows),
		_cols(cols),
		_data(rows*cols, 0)
	{
		if (!rows || !cols)
		{
			throw "invalid matrix size\n";
		}
	}


	template<typename T>
	unsigned long matrix<T>::rows() const
	{
		return _rows;
	}


	template<typename T>
	unsigned long matrix<T>::cols() const
	{
		return _cols;
	}


	template<typename T>
	std::vector<T> matrix<T>::data() const
	{
		return _data;
	}

	
	template<typename T>
	matrix<T> matrix<T>::identity(unsigned long size)
	{
		matrix<T> result(size, size);

		long count = 0;
		std::generate(result._data.begin(), result._data.end(),
			[&count, size]()
		{
			return !(count++ % (size + 1));
		});

		return result;
	}


	template<typename T>
	T& matrix<T>::operator()(unsigned long row, unsigned long col) const
	{
		if (row >= _rows || col >= _cols)
		{
			throw "position out of range\n";
		}

		return (T&)_data[col + _cols*row];
	}


	template<typename T>
	matrix<T> matrix<T>::operator*(const matrix<T>& other) const
	{
		if (_cols != other._rows)
		{
			throw "matrix dimensions are not multiplicable\n";
		}

		matrix<T> result(_rows, other._cols);

		for (unsigned long r = 0; r < _rows; r++)
		{
			for (unsigned long ocol = 0; ocol < other._cols; ocol++)
			{
				for (unsigned long c = 0; c < _cols; c++)
				{
					result(r, ocol) += (*this)(r, c) * other(c, ocol);
				}
			}
		}

		return result;
	}


	template<typename T>
	matrix<T> matrix<T>::transpose() const
	{
		matrix<T> result(_cols, _rows);

		for (unsigned long r = 0; r < _rows; ++r)
		{
			for (unsigned long c = 0; c < _cols; ++c)
			{
				result(c, r) += (*this)(r, c);
			}
		}
		return result;
	}
}

#endif