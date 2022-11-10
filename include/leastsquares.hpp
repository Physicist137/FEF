#pragma once
#include <Eigen/Dense>
#include <vector>
#include <functional>

template <typename T>
class LeastSquaresData {
	T a;
	T b;
	T aberr;
	T reerr;

public:
	LeastSquaresData(const T& lin, const T& ind, const T& ab, const T& re) 
	: a(lin), b(ind), aberr(ab), reerr(re) {}
	
	const T& independent() const {return b;}
	const T& linear() const {return a;}
	const T& absolute_quadratic_error() const {return aberr;}
	const T& relative_quadratic_error() const {return reerr;}
};


template <typename T>
T identity(const T& value) {
	return value;
}

template <typename T>
LeastSquaresData<T> leastSquares(const std::vector<T>& x, const std::vector<T>& y,
const std::function<T(T)>& fy = identity<T>, const std::function<T(T)>& fx = identity<T>) {
	unsigned size = x.size();
	T tsize = static_cast<T>(size);

	T avx = 0.0;
	T avy = 0.0;
	for (unsigned i = 0; i < size; ++i) {
		avx += fx(x[i]) / tsize;
		avy += fy(y[i]) / tsize;		
	}

	T num = 0.0;
	T den = 0.0;
	for (unsigned i = 0; i < size; ++i) {
		num += fx(x[i]) * (fy(y[i]) - avy);
		den += fx(x[i]) * (fx(x[i]) - avx);
	}
	
	T linear = num / den;
	T independent = avy - linear * avx;

	T abs_error = 0.0;
	T rel_error = 0.0;
	for (unsigned i = 0; i < size; ++i) {
		T predicted = independent + linear * x[i];
		T diff = y[i] - predicted;

		abs_error += diff * diff;
		rel_error += diff * diff / (y[i] * y[i]);
	}

	return LeastSquaresData<T>(linear, independent, abs_error, rel_error);
}

template <typename T>
LeastSquaresData<T> leastSquares(
const Eigen::Matrix<T, Eigen::Dynamic, 1>& x, const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
const std::function<T(T)>& fy = identity<T>, const std::function<T(T)>& fx = identity<T>) {
	unsigned size = x.size();
	T tsize = static_cast<T>(size);

	T avx = 0.0;
	T avy = 0.0;
	for (unsigned i = 0; i < size; ++i) {
		avx += fx(x(i)) / tsize;
		avy += fy(y(i)) / tsize;		
	}

	T num = 0.0;
	T den = 0.0;
	for (unsigned i = 0; i < size; ++i) {
		num += fx(x(i)) * (fy(y(i)) - avy);
		den += fx(x(i)) * (fx(x(i)) - avx);
	}
	
	T linear = num / den;
	T independent = avy - linear * avx;

	T abs_error = 0.0;
	T rel_error = 0.0;
	for (unsigned i = 0; i < size; ++i) {
		T predicted = independent + linear * x(i);
		T diff = y(i) - predicted;

		abs_error += diff * diff;
		rel_error += diff * diff / (y(i) * y(i));
	}

	return LeastSquaresData<T>(linear, independent, abs_error, rel_error);
}

