//
// Copyright (c) 2003 Joao Abecasis
//

#ifndef _agedo_LINEAR_ALGEBRA_H
#define _agedo_LINEAR_ALGEBRA_H

#include <cmath>

namespace LSODA { namespace LinearAlgebra {

	// absolute value
	template<typename value_type> inline value_type Absolute(const value_type &val) { return val < value_type(0) ? -val : val; }
	template<> inline int Absolute<int>(const int &val) { return std::abs(val); }
//	template<> inline float Absolute<float>(const float &val) { return std::fabsf(val); }
	template<> inline double Absolute<double>(const double &val) { return std::fabs(val); }

	// weighted max-norm of vector of length n
	//   = max(i=0...n-1)( abs(vector[i]) * weight[i] )
	template<typename value_type, typename size_type>
	inline value_type VectorMaxNorm(size_type n, const value_type * const vector, const value_type * const weight)
	{
		--n;
		value_type temp, vmax = Absolute<value_type>(vector[n])*weight[n];

		while (--n >= 0) {
			temp = Absolute(vector[n]) * weight[n];
			temp > vmax ? vmax = temp : 0;
		}
		return vmax;
	}

	// weighted max-norm of n by n matrix
	//   = max(i=0...n-1)( weight[i] * sum(j=0...n-1)( abs(matrix[i + j*n]) / weight[j] ) )
	template<typename value_type, typename size_type>
	inline value_type MatrixMaxNorm(const size_type n, const value_type * matrix, const value_type * const weight)
	{
		size_type i = 0, j;
		value_type temp, mmax = 0;

		for (; i < n; ++i, ++matrix) {
			for (j = 0, temp = value_type(0); j < n ; ++j) {
				temp += Absolute(matrix[j*n]) / weight[j];
			}
			temp *= weight[i];
			temp > mmax ? mmax = temp : 0;
		}
		return mmax;
	}

	// max-norm of banded n by n matrix
	// This function computes the norm of a banded N by N matrix,
	// stored in the array MATRIX, that is consistent with the weighted max-norm
	// on vectors, with weights stored in the array WEIGHT.
	// ML and MU are the lower and upper half-bandwidths of the matrix.
	// NRA is the first dimension of the A array, NRA >= ML+MU+1.
	// In terms of the matrix elements a[i,j], the norm is given by:
	//   = max(i=0...n-1) ( weight[i] * sum(j=0...n-1)( abs( matrix[i,j] ) / weight[j] ) )
	template<typename value_type, typename size_type>
	inline value_type BandedMatrixMaxNorm(const size_type n, const value_type * const matrix, const size_type nRows, const size_type ml, const size_type mu, const value_type * const weight)
	{
		size_type i = 0, j, jHi;
		value_type temp, mmax = 0;

		for (; i < n; ++i) {
			temp = value_type(0);
			jHi = std::min(i + mu + 1, n);
			for (j = std::max(i - ml, size_type(0)); j < jHi; ++j) {
				temp += Absolute(matrix[i + mu - j + j * nRows]) / weight[j];
			}
			temp *= weight[i];
			temp > mmax ? mmax = temp: 0 ;
		}
		return mmax;
	}

	// find the smallest (1-based) index of that component of a vector having the maximum magnitude
	//    = min(i=0...n-1) : vector[i] == max(j=0...n-1)( abs(vector[j-1]) )
	//    returns 0 if n <= 0;
	template<typename value_type, typename size_type>
	inline size_type IndexOfMax(const size_type n, const value_type * vector)
	{
		size_type ret_val = 0;
		if (n <= 0) return ret_val;

		ret_val = 1;
		if (n == 1) return ret_val;

		--vector;
		value_type temp, vmax = Absolute(vector[1]);

		for (size_type i = 1; i <= n; ++i)
		{
			temp = Absolute(vector[i]);
			temp > vmax ? (ret_val = i, vmax = temp) : 0;
		}

		return ret_val;
	}

	// compute a constant times a vector
	//   sets : (i=0...n-1) ( X[i] = A * X[i] )
	template<typename value_type, typename size_type>
	inline void AX(size_type n, const value_type &A, value_type * const X)
	{
		for (; n > 3; n -= 4)
		{
			X[n-1] *= A;
			X[n-2] *= A;
			X[n-3] *= A;
			X[n-4] *= A;
		}

		while (--n >= 0)
		{
			X[n] *= A;
		}

		return;
	}

	// compute a constant times a vector plus a vector
	//   sets : (i=0...n-1)( Y[i] = A * X[i] + Y[i] )
	template<typename value_type, typename size_type>
	inline void AXplusY(size_type n, const value_type &A, const value_type * const X, value_type * const Y)
	{
		for (; n > 3; n -= 4)
		{
			Y[n-1] += A * X[n-1];
			Y[n-2] += A * X[n-2];
			Y[n-3] += A * X[n-3];
			Y[n-4] += A * X[n-4];
		}

		while (--n >= 0)
		{
			Y[n] += A * X[n];
		}

		return;
	}

	// compute the inner product of two vectors
	//   = sum(i=1...n)( X[i] * Y[i] )
	template<typename value_type, typename size_type>
	inline value_type DotProduct(size_type n, const value_type * const X, const value_type * const Y)
	{
		value_type ret_val = 0;

		for (; n > 3; n -= 4)
		{
			ret_val += X[n-1] * Y[n-1] + X[n-2] * Y[n-2] + X[n-3] * Y[n-3] + X[n-4] * Y[n-4];
		}

		while (--n >= 0)
		{
			ret_val += X[n] * Y[n];
		}

		return ret_val;
	}
} /* LinearAlgebra */ } /* LSODA */

#endif
