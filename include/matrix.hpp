#pragma once
#include <math/vector.hpp>
#include <array>
#include <cmath>
#include <initializer_list>
#include <algorithm>
#include <iostream>

/*
To think off:
Make: (i, j) == [i][j].
Same with the const values.
*/

namespace math {

/// Echelon Reduction type enumeration
enum class ReductionType {
	UpperLeft, UpperRight,
	LowerLeft, LowerRight
};

template <typename T, unsigned Rows, unsigned Cols>
class Matrix {
    // Matrix formed by column-vectors:
    // std::array< math::vector<T, Rows>, Cols> _data; Is equivalent:
    std::array< math::Vector<T, Rows>, Cols> _data;
    static_assert(Cols*Rows > 0, "Can't make matrix with no cells");

public:
    Matrix();
    Matrix(std::initializer_list<math::Vector<T, Cols>> list);
	Matrix(std::initializer_list<T> list);

    Matrix(const Matrix<T, Rows, Cols>& other);
    Matrix(Matrix<T, Rows, Cols>& other);

    const T& operator()(unsigned i, unsigned j) const;
    T& operator()(unsigned i, unsigned j);

    // This must be changed:
    /// Returns the jth col-vector of the Matrix. Notice. j < Cols.
    const math::Vector<T, Rows>& operator[](unsigned j) const;
    math::Vector<T, Rows>& operator[](unsigned j);

    Matrix operator+() const;
    Matrix operator-() const;

    Matrix& operator+=(const Matrix& mat);
    Matrix& operator-=(const Matrix& mat);
    Matrix& operator*=(const Matrix& mat);

    Matrix& operator*=(const T& value);
    Matrix& operator/=(const T& value);

    static const Matrix<T, Rows, Cols> zeros();
    static const Matrix<T, Rows, Cols> ones();
    static const Matrix<T, Rows, Cols> eye();

    Vector<T, Cols> row(unsigned i) const;
    Vector<T, Rows> col(unsigned j) const;

    Matrix<T, Cols, Rows> transpost() const;

    // Reduced Row Echelon Form
    Matrix<T, Rows, Cols> rref(ReductionType type) const;
    Matrix<T, Rows, Cols> rref(T& value, ReductionType type) const;
    Matrix<T, Rows, Cols> rref(Vector<T, Rows>& value, ReductionType type) const;
    Matrix<T, Rows, Cols> rref(Matrix<T, Rows, Cols>& value, ReductionType type) const;

    template <typename ReductionHelper>
    Matrix<T, Rows, Cols> rref(ReductionHelper&& helper, ReductionType type) const;

	// TODO: Decide!
	// 	Shall we specialize Matrix class for square ones and put det and inverse?
	//	Shall we put det/inverse here?
	T det() const;
	Matrix<T, Rows, Rows> inverse() const;

    // Other Functions. Swap, etc.
	void swapline(unsigned a, unsigned b);
};

// Definition of alias -----------------------------------------------------------
using Matrix2 = Matrix<double, 2, 2>;
using Matrix3 = Matrix<double, 3, 3>;
using Matrix4 = Matrix<double, 4, 4>;


// Declaration of functions: -----------------------------------------------------
template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator+(const Matrix<T, Rows, Cols>& matrix, const Matrix<T, Rows, Cols>& other);

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator-(const Matrix<T, Rows, Cols>& matrix, const Matrix<T, Rows, Cols>& other);

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator*(const Matrix<T, Rows, Cols>& matrix, const T& value);

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator*(const T& value, const Matrix<T, Rows, Cols>& matrix);

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator/(const Matrix<T, Rows, Cols>& matrix, const T& value);

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator/(const T& value, const Matrix<T, Rows, Cols>& matrix);

template <typename T, unsigned Rows, unsigned Inter, unsigned Cols>
Matrix<T, Rows, Cols> operator*(const Matrix<T, Rows, Inter>& matrix, const Matrix<T, Inter, Cols>& other);

template <typename T, unsigned Rows, unsigned Cols>
Vector<T, Rows> operator*(const Matrix<T, Rows, Cols>& matrix, const Vector<T, Cols>& vec);

// Definition of member-functions: -----------------------------------------------
template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols>::Matrix() {
    // All vectors will be constructed and initialized to T()
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols>::Matrix(std::initializer_list<math::Vector<T, Cols>> list) {
    for (unsigned j = 0; j < Cols; ++j) {
        std::copy(list.begin(), list.end(), _data.begin());
    }
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols>::Matrix(std::initializer_list<T> list) {
	unsigned i = 0;
	for (const T& value : list) {
		this->operator()(i / Cols, i % Cols) = value;
		++i;
	}
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols>::Matrix(const Matrix<T, Rows, Cols>& other) : _data(other._data) {

}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols>::Matrix(Matrix<T, Rows, Cols>& other) : _data(other._data) {

}

template <typename T, unsigned Rows, unsigned Cols>
const T& Matrix<T, Rows, Cols>::operator()(unsigned i, unsigned j) const {
    #ifdef DEBUG
        if (i >= Rows  || j >= Cols) throw std::logic_error("Invalid index")
    #endif

    return _data[j][i];
}

template <typename T, unsigned Rows, unsigned Cols>
T& Matrix<T, Rows, Cols>::operator()(unsigned i, unsigned j) {
    #ifdef DEBUG
        if (i >= Rows  || j >= Cols) throw std::logic_error("Invalid index")
    #endif

    return _data[j][i];
}

template <typename T, unsigned Rows, unsigned Cols>
const math::Vector<T, Rows>& Matrix<T, Rows, Cols>::operator[](unsigned j) const {
    #ifdef DEBUG
        if (j >= Cols) throw std::logic_error("Invalid index")
    #endif

    return _data[j];
}

template <typename T, unsigned Rows, unsigned Cols>
math::Vector<T, Rows>& Matrix<T, Rows, Cols>::operator[](unsigned j) {
    #ifdef DEBUG
        if (j >= Cols) throw std::logic_error("Invalid index")
    #endif

    return _data[j];
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator+() const {
    return *this;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::operator-() const {
    Matrix<T, Rows, Cols> result;
    for (unsigned i = 0; i < Cols; ++i) result[i] = -_data[i];
    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols>& Matrix<T, Rows, Cols>::operator+=(const Matrix& mat) {
    for (unsigned i = 0; i < Cols; ++i) _data[i] += mat[i];
    return *this;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols>& Matrix<T, Rows, Cols>::operator-=(const Matrix& mat) {
    for (unsigned i = 0; i < Cols; ++i) _data[i] -= mat[i];
    return *this;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols>& Matrix<T, Rows, Cols>::operator*=(const Matrix& mat) {
    // TODO.
    return (*this) = (*this) * mat;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols>& Matrix<T, Rows, Cols>::operator*=(const T& value) {
    for (unsigned i = 0; i < Cols; ++i) _data[i] *= value;
    return *this;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols>& Matrix<T, Rows, Cols>::operator/=(const T& value) {
    for (unsigned i = 0; i < Cols; ++i) _data[i] /= value;
    return *this;
}

template <typename T, unsigned Rows, unsigned Cols>
const Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::zeros() {
    return Matrix<T, Rows, Cols>();
}

template <typename T, unsigned Rows, unsigned Cols>
const Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::ones() {
    Matrix<T, Rows, Cols> result;
    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            result._data[j][i] = 1;
        }
    }

    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
const Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::eye() {
    Matrix<T, Rows, Cols> result;
    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            result._data[j][i] = (i == j);
        }
    }

    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
Vector<T, Cols> Matrix<T, Rows, Cols>::row(unsigned i) const {
    Vector<T, Cols> result;
    for (unsigned j = 0; j < Cols; ++j) result(j) = this->operator()(i, j);
    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
Vector<T, Rows> Matrix<T, Rows, Cols>::col(unsigned j) const {
    Vector<T, Rows> result;
    for (unsigned i = 0; j < Rows; ++i) result(i) = this->operator()(i, j);
    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Cols, Rows> Matrix<T, Rows, Cols>::transpost() const {
    Matrix<T, Cols, Rows> result;
    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            result(j, i) = this->operator()(i, j);
        }
    }

    return result;
}

// Definition of non-member-functions: -------------------------------------------
template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator+(const Matrix<T, Rows, Cols>& matrix, const Matrix<T, Rows, Cols>& other) {
    Matrix<T, Rows, Cols> result;
    for (unsigned i = 0; i < Cols; ++i) result[i] = matrix[i] + other[i];
    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator-(const Matrix<T, Rows, Cols>& matrix, const Matrix<T, Rows, Cols>& other) {
    Matrix<T, Rows, Cols> result;
    for (unsigned i = 0; i < Cols; ++i) result[i] = matrix[i] - other[i];
    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator*(const Matrix<T, Rows, Cols>& matrix, const T& value) {
    Matrix<T, Rows, Cols> result;
    for (unsigned i = 0; i < Cols; ++i) result[i] = matrix[i] * value;
    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator*(const T& value, const Matrix<T, Rows, Cols>& matrix) {
    return matrix * value;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator/(const Matrix<T, Rows, Cols>& matrix, const T& value) {
    Matrix<T, Rows, Cols> result;
    for (unsigned i = 0; i < Cols; ++i) result[i] = matrix[i] / value;
    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> operator/(const T& value, const Matrix<T, Rows, Cols>& matrix) {
    return matrix / value;
}

template <typename T, unsigned Rows, unsigned Inter, unsigned Cols>
Matrix<T, Rows, Cols> operator*(const Matrix<T, Rows, Inter>& matrix, const Matrix<T, Inter, Cols>& other) {
    Matrix<T, Rows, Cols> result;
    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            for (unsigned k = 0; k < Inter; ++k) {
                result(i, j) += matrix(i, k) * other(k, j);
            }
        }
    }

    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
Vector<T, Rows> operator*(const Matrix<T, Rows, Cols>& matrix, const Vector<T, Cols>& vec) {
    Vector<T, Rows> result;

    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            result(i) += matrix(i, j) * vec(j);
        }
    }

    return result;
}

template <typename T, unsigned Rows, unsigned Cols>
void Matrix<T, Rows, Cols>::swapline(unsigned a, unsigned b) {
	// TODO: Fix me for direct accessment to _data, as soon as _data[i][j] is working.
	for (unsigned j = 0; j < Cols; ++j) std::swap(this->operator()(a, j), this->operator()(b, j));
}

// Internal namespace for reduction helpers
namespace internal {

template <typename T, unsigned Rows, unsigned Cols>
struct ReductionHelperWithMatrix {
	Matrix<T, Rows, Cols>& mat;

	ReductionHelperWithMatrix(Matrix<T, Rows, Cols>& s) : mat(s) {}
	inline void applyLineSwap(unsigned a, unsigned b) {
		mat.swapline(a, b);
		//for (unsigned j = 0; j < Cols; ++j) std::swap(mat(a, j), mat(b, j));
	}

	inline void applySubtractLines(unsigned a, unsigned b) {
		for (unsigned j = 0; j < Cols; ++j) mat(a, j) -= mat(b, j);
	}

	inline void applyScalar(unsigned line, const T& scalar) {
		for (unsigned j = 0; j < Cols; ++j) mat(line, j) *= scalar;
	}
};


template <typename T, unsigned D>
struct ReductionHelperWithVector {
	Vector<T, D>& vec;

	ReductionHelperWithVector(Vector<T, D>& s) : vec(s) {}
	inline void applyLineSwap(unsigned a, unsigned b) {std::swap(vec[a], vec[b]);}
	inline void applySubtractLines(unsigned a, unsigned b) {vec[a] -= vec[b];}
	inline void applyScalar(unsigned line, const T& scalar) {vec[line] *= scalar;}
};


template <typename T>
struct ReductionHelperWithScalar {
	T& scalar;

	ReductionHelperWithScalar(T& s) : scalar(s) {}
	inline void applyLineSwap(unsigned, unsigned) {scalar *= static_cast<T>(-1);}
	inline void applySubtractLines(unsigned, unsigned) {}
	inline void applyScalar(unsigned, const T& s) {scalar *= s;}
};

template <typename T>
struct ReductionHelperWithNothing {
	inline void applyLineSwap(unsigned, unsigned) {}
	inline void applySubtractLines(unsigned, unsigned) {}
	inline void applyScalar(unsigned, const T&) {}
};
} // internal namespace

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::rref(ReductionType type) const {
	return rref(internal::ReductionHelperWithNothing<T>(), type);
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::rref(T& value, ReductionType type) const {
	return rref(internal::ReductionHelperWithScalar<T>(value), type);
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::rref(Vector<T, Rows>& value, ReductionType type) const {
	return rref(internal::ReductionHelperWithVector<T, Rows>(value), type);
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::rref(Matrix<T, Rows, Cols>& value, ReductionType type) const {
	return rref(internal::ReductionHelperWithMatrix<T, Rows, Cols>(value), type);
}


template <typename T, unsigned Rows, unsigned Cols>
template <typename ReductionHelper>
Matrix<T, Rows, Cols> Matrix<T, Rows, Cols>::rref(ReductionHelper&& helper, ReductionType type) const {
	Matrix<T, Rows, Cols> result(*this);
	internal::ReductionHelperWithMatrix<T, Rows, Cols> mhelper(result);

	for (unsigned k = 0; k < Rows-1  or  k < Cols-1; ++k) {
		// Null pivot indicator with no chance of swap
		bool nullpivot = false;

		// Pivot matrix indexes. Row, Col.
		unsigned ip, jp;

		// Selecting pivot indexes.
		switch (type) {
			case ReductionType::LowerRight: ip=Cols-1-k;  jp=k;         break;	// Pivot on Right-superior corner
			case ReductionType::LowerLeft:  ip=k;         jp=k;         break;  // Pivot on Left-superior corner
			case ReductionType::UpperLeft:  ip=k;         jp=Rows-1-k;  break;	// Pivot on Left-inferior corner
			case ReductionType::UpperRight: ip=Cols-1-k;  jp=Rows-1-k;	break;	// Pivot on Right-inferior corner
		}

		// The pivot value
		T pivot = result(ip, jp);

		// If pivot is zero, iterate until find non zero, and swap.
		if (pivot == 0) {
			nullpivot = true;
			for (unsigned i = k+1; i < Rows; ++i) {
				// Cell matrix indexes. Row, Col.
				unsigned ic = ip;
				unsigned jc = jp;

				// Selecting cell indexes.
				switch (type) {
					case ReductionType::LowerRight: ic = i; break;
					case ReductionType::LowerLeft:  ic = i; break;
					case ReductionType::UpperLeft:  ic = Rows-i-1; break;
					case ReductionType::UpperRight: ic = Rows-i-1; break;
				}

				// The cell value
				T cell = result(ic, jc);

				// If cell is nonzero, swap.
				if (cell != 0) {
					helper.applyLineSwap(ip, ic);
					mhelper.applyLineSwap(ip, ic);
					nullpivot = false;
					break;
				}
			}
		}

		// If entire pivot line is null, move on...
		if (nullpivot) continue;
		// FIXME. (just increase the pivot-j without passing to other line).
		// But then, it won't work for determinants, and inverse, in case this happen... =(

		// Reduce the matrix
		for (unsigned i = k+1; i < Rows; ++i) {
			// Cell matrix indexes. Row, Col.
			unsigned ic = ip;
			unsigned jc = jp;

			// Select the cell indexes.
			switch (type) {
				case ReductionType::LowerRight: ic = i; break;
				case ReductionType::LowerLeft:  ic = i; break;
				case ReductionType::UpperLeft:  ic = Rows-i-1; break;
				case ReductionType::UpperRight: ic = Rows-i-1; break;
			}

			// Cell value
			T cell = result(ic, jc);

			// Already reduced. Move on...
			if (cell == 0) continue;

			// Apply reduction
			T value = pivot / cell;
			helper.applyScalar(ic, value);
			helper.applySubtractLines(ic, ip);
			mhelper.applyScalar(ic, value);
			mhelper.applySubtractLines(ic, ip);
		}
	}

	return result;
}

template <typename T, unsigned Rows, unsigned Cols>
T Matrix<T, Rows, Cols>::det() const {
	static_assert(Rows == Cols, "Matrix must be squared for det() to work");

	T result = 1;
	Matrix reduced = this->rref(result, ReductionType::LowerLeft);

	T diag = 1;
	for (unsigned i = 0; i < Rows; ++i) diag *= reduced(i, i);
	return diag / result;
}

template <typename T, unsigned Rows, unsigned Cols>
Matrix<T, Rows, Rows> Matrix<T, Rows, Cols>::inverse() const {
	static_assert(Rows == Cols, "Matrix must be squared for det() to work");

	Matrix identity = Matrix<T, Rows, Cols>::eye();
	Matrix firstReduced = this->rref(identity, ReductionType::LowerLeft);
	Matrix secondReduced = firstReduced.rref(identity, ReductionType::UpperRight);

	double d = det();
	if (d == 0) throw std::logic_error("This matrix is non-invertible");

	for (unsigned i = 0; i < Rows; ++i) {
		for (unsigned j = 0; j < Cols; ++j) identity(i, j) /= secondReduced(i, i);
	}

	return identity;
}

/*
for (unsigned i = 0; i < Rows; ++i) {
	for (unsigned j = 0; j < Cols; ++j) std::cout << this->operator()(i, j) << " ";
	std::cout << std::endl;
}
std::cout << std::endl;
*/


} // Math namespace
