#ifndef _APPROX_HPP_
#define _APPROX_HPP_

#include <vector>
#include "matrix.hpp"
#include "givensQR.hpp"

namespace nma
{
	template<typename T>
	std::vector<T> polyfit(const std::vector<T>& X, const std::vector<T>& Y, long power)
	{
		if (X.size() != Y.size())
			throw "X and Y vector sizes do not match\n";

		power++;

		size_t count = X.size();
		matrix<T> XMatrix(count, power);
		matrix<T> YMatrix(count, 1);

		for (size_t i = 0; i < count; i++)
		{
			YMatrix(i, 0) = Y[i];
		}

		for (size_t row = 0; row < count; row++)
		{
			T nVal = 1.0f;
			for (long col = 0; col < power; col++)
			{
				XMatrix(row, col) = nVal;
				nVal *= X[row];
			}
		}

		// transpose X matrix
		matrix<T> XtMatrix(XMatrix.transpose());
		// multiply transposed X matrix with X matrix
		matrix<T> XtXMatrix(XtMatrix * XMatrix);
		// multiply transposed X matrix with Y matrix
		matrix<T> XtYMatrix(XtMatrix * YMatrix);

		Givens<T> oGivens;
		oGivens.decompose(XtXMatrix);
		matrix<T> coefs = oGivens.solve(XtYMatrix);

		return coefs.data();
	}


	template<typename T>
	std::vector<T> polyval(const std::vector<T>& coefs, const std::vector<T>& X)
	{
		size_t count = X.size();
		size_t power = coefs.size();
		std::vector<T> Y(count);

		for (size_t i = 0; i < count; i++)
		{
			T nY = 0;
			T nXT = 1;
			T nX = X[i];
			for (size_t j = 0; j < power; j++)
			{
				nY += coefs[j] * nXT;
				nXT *= nX;
			}
			Y[i] = nY;
		}

		return Y;
	}

}

#endif