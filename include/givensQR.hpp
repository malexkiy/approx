#ifndef _GIVENSQR_HPP_
#define _GIVENSQR_HPP_

#include "matrix.hpp"
#include <cmath>
#include <algorithm>

namespace nma
{
	template<typename T>
	class Givens
	{
	private:
		matrix<T> _oQ, _oR, _oJ;

		void givensRotation(T a, T b);
		void preMultiplyGivens(matrix<T>& mtx, long i, long j) const;

	public:
		Givens() : _oJ(2, 2), _oQ(1, 1), _oR(1, 1) {}
		
		const matrix<T> inverse(matrix<T>& mtx) const;

		void decompose(matrix<T>& mtx);
		matrix<T> solve(matrix<T>& mtx) const;

		const matrix<T>& GetQ() const;
		const matrix<T>& GetR() const;
	};


	template<typename T>
	void Givens<T>::givensRotation(T a, T b)
	{
		T t, s, c;

		if (b == 0)
		{
			c = (a >= 0) ? 1 : -1;
			s = 0;
		}
		else if (a == 0)
		{
			c = 0;
			s = (b >= 0) ? -1 : 1;
		}
		else if (abs(b) > abs(a))
		{
			t = a / b;
			s = -1 / sqrt(1 + t*t);
			c = -s*t;
		}
		else
		{
			t = b / a;
			c = 1 / sqrt(1 + t*t);
			s = -c*t;
		}
		_oJ(0, 0) = c;
		_oJ(0, 1) = -s;

		_oJ(1, 0) = s;
		_oJ(1, 1) = c;
	}


	template<typename T>
	void Givens<T>::preMultiplyGivens(matrix<T>& mtx, long i, long j) const
	{
		long rowSize = mtx.cols();

		for (long row = 0; row < rowSize; row++)
		{
			double temp = mtx(i, row) * _oJ(0, 0) + mtx(j, row) * _oJ(0, 1);
			mtx(j, row) = mtx(i, row) * _oJ(1, 0) + mtx(j, row) * _oJ(1, 1);
			mtx(i, row) = temp;
		}
	}


	template<typename T>
	const matrix<T> Givens<T>::inverse(matrix<T>& mtx) const
	{
		if (mtx.cols() != mtx.rows())
		{
			throw "matrix has to be square\n";
		}

		matrix<T> identity = matrix<T>::identity(mtx.rows());
		decompose(mtx);

		return solve(identity);
	}


	template<typename T>
	void Givens<T>::decompose(matrix<T>& mtx)
	{
		long rows = mtx.rows();
		long cols = mtx.cols();


		if (rows == cols)
		{
			cols--;
		}
		else if (rows < cols)
		{
			cols = rows - 1;
		}

		_oQ = matrix<T>::identity(rows);
		_oR = mtx;

		for (long j = 0; j < cols; j++)
		{
			for (long i = j + 1; i < rows; i++)
			{
				givensRotation(_oR(j, j), _oR(i, j));
				preMultiplyGivens(_oR, j, i);
				preMultiplyGivens(_oQ, j, i);
			}
		}

		_oQ = _oQ.transpose();
	}


	template<typename T>
	matrix<T> Givens<T>::solve(matrix<T>& mtx) const
	{
		matrix<T> oQtM(_oQ.transpose() * mtx);
		long cols = _oR.cols();
		matrix<T> oS(1, cols);

		for (long i = cols - 1; i >= 0; i--)
		{
			oS(0, i) = oQtM(i, 0);
			for (long j = i + 1; j < cols; j++)
			{
				oS(0, i) -= oS(0, j) * _oR(i, j);
			}
			oS(0, i) /= _oR(i, i);
		}

		return oS;
	}


	template<typename T>
	const matrix<T>& Givens<T>::GetQ() const
	{
		return _oQ;
	}


	template<typename T>
	const matrix<T>& Givens<T>::GetR() const
	{
		return _oR;
	}
}

#endif