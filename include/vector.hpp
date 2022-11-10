#pragma once
#include <array>
#include <stdexcept>
#include <cmath>
#include <initializer_list>
#include <algorithm>
#include <ostream>

namespace math {
template <typename T, unsigned D>
class Vector {
    std::array<T, D> _data;
    static_assert(D > 0, "Can't make a vector with no dimentions");

public:
    Vector();
    Vector(std::initializer_list<T> list);

    Vector(const Vector<T, D>& other);
    Vector(Vector<T, D>& other);

    const T& operator[](unsigned i) const;
    T& operator[](unsigned i);

    const T& operator()(unsigned i) const;
    T& operator()(unsigned i);

    T x() const;
    T y() const;
    T z() const;
    T t() const;

    Vector operator+() const;
    Vector operator-() const;
    Vector operator/(const T& value);

    Vector& operator+=(const Vector& vec);
    Vector& operator-=(const Vector& vec);
    Vector& operator*=(const Vector& vec);
    Vector& operator/=(const Vector& vec);

    Vector& operator*=(const T& value);
    Vector& operator/=(const T& value);

	Vector cross(const Vector& value) const;
    Vector unit() const;

    T dot(const Vector& vec) const;
    T dot() const;

    template <typename E=T>
    E length() const;
};

// Definition of alias -----------------------------------------------------------
using Vector2 = Vector<double, 2>;
using Vector3 = Vector<double, 3>;
using Vector4 = Vector<double, 4>;

// Declaration of functions: -----------------------------------------------------
template <typename T, unsigned D> Vector<T, D> operator+(const Vector<T, D>& vec, const Vector<T, D>& other);
template <typename T, unsigned D> Vector<T, D> operator-(const Vector<T, D>& vec, const Vector<T, D>& other);
template <typename T, unsigned D> Vector<T, D> operator*(const Vector<T, D>& vec, const Vector<T, D>& other);
template <typename T, unsigned D> Vector<T, D> operator/(const Vector<T, D>& vec, const Vector<T, D>& other);

template <typename T, unsigned D> Vector<T, D> operator*(const Vector<T, D>& vec, const T& value);
template <typename T, unsigned D> Vector<T, D> operator*(const T& value, const Vector<T, D>& vec);
template <typename T, unsigned D> Vector<T, D> operator/(const Vector<T, D>& vec, const T& value);
template <typename T, unsigned D> Vector<T, D> operator/(const T& value, const Vector<T, D>& vec);

template <typename T, unsigned D> std::ostream& operator<<(std::ostream& os, const math::Vector<T, D>& vec);

// Definition of memberfunctions: ------------------------------------------------
template <typename T, unsigned D>
inline Vector<T, D>::Vector() {
    for (unsigned i = 0; i < D; ++i) _data[i] = T(0);
}

template <typename T, unsigned D>
inline Vector<T, D>::Vector(std::initializer_list<T> list) {
    std::copy(list.begin(), list.end(), _data.begin());
}

template <typename T, unsigned D>
inline Vector<T, D>::Vector(const Vector<T, D>& other) : _data(other._data) {

}

template <typename T, unsigned D>
inline Vector<T, D>::Vector(Vector<T, D>& other) : _data(other._data) {

}

template <typename T, unsigned D>
T Vector<T, D>::x() const {
    static_assert(D>=1, "The x-component of this vector does not exist");
    return _data[0];
}

template <typename T, unsigned D>
T Vector<T, D>::y() const {
    static_assert(D>=2, "The y-component of this vector does not exist");
    return _data[1];
}

template <typename T, unsigned D>
T Vector<T, D>::z() const {
    static_assert(D>=3, "The z-component of this vector does not exist");
    return _data[2];
}

template <typename T, unsigned D>
T Vector<T, D>::t() const {
    static_assert(D>=4, "The t-component of this vector does not exist");
    return _data[3];
}

template <typename T, unsigned D>
inline const T& Vector<T, D>::operator[](unsigned i) const {
    #ifdef DEBUG
        if (i >= D) throw std::logic_error("Invalid index")
    #endif

    return _data[i];
}

template <typename T, unsigned D>
inline T& Vector<T, D>::operator[](unsigned i) {
    #ifdef DEBUG
        if (i >= D) throw std::logic_error("Invalid index")
    #endif

    return _data[i];
}

template <typename T, unsigned D>
inline const T& Vector<T, D>::operator()(unsigned i) const {
    #ifdef DEBUG
        if (i >= D) throw std::logic_error("Invalid index")
    #endif
    return this->operator[](i);
}

template <typename T, unsigned D>
inline T& Vector<T, D>::operator()(unsigned i) {
    #ifdef DEBUG
        if (i >= D) throw std::logic_error("Invalid index")
    #endif
    return this->operator[](i);
}

template <typename T, unsigned D>
inline Vector<T, D> Vector<T, D>::operator+() const {
    return *this;
}

template <typename T, unsigned D>
inline Vector<T, D> Vector<T, D>::operator-() const {
    Vector vec;
    for (unsigned i = 0; i < D; ++i) vec[i] = -_data[i];
    return vec;
}

template <typename T, unsigned D>
inline Vector<T, D> Vector<T, D>::operator/(const T& value) {
    #ifdef DEBUG
    	for (unsigned i = 0; i < D; ++i)
    		if (_data[i] == 0)
    			throw std::logic_error("Can't divide by zero");
    #endif

    Vector vec;
    for (unsigned i = 0; i < D; ++i) vec[i] = _data[i] / value;
    return vec;
}

template <typename T, unsigned D>
inline Vector<T, D>& Vector<T, D>::operator+=(const Vector<T, D>& vec) {
    for (unsigned i = 0; i < D; ++i) _data[i] += vec[i];
    return *this;
}

template <typename T, unsigned D>
inline Vector<T, D>& Vector<T, D>::operator-=(const Vector<T, D>& vec) {
    for (unsigned i = 0; i < D; ++i) _data[i] -= vec[i];
    return *this;
}

template <typename T, unsigned D>
inline Vector<T, D>& Vector<T, D>::operator*=(const Vector<T, D>& vec) {
    for (unsigned i = 0; i < D; ++i) _data[i] *= vec[i];
    return *this;
}

template <typename T, unsigned D>
inline Vector<T, D>& Vector<T, D>::operator/=(const Vector<T, D>& vec) {
    #ifdef DEBUG
    	for (unsigned i = 0; i < D; ++i)
    		if (vec._data[i] == 0)
    			throw std::logic_error("Can't divide by zero");
    #endif

    for (unsigned i = 0; i < D; ++i) _data[i] /= vec[i];
    return *this;
}

template <typename T, unsigned D>
inline Vector<T, D>& Vector<T, D>::operator*=(const T& value) {
    for (unsigned i = 0; i < D; ++i) _data[i] *= value;
    return *this;
}

template <typename T, unsigned D>
inline Vector<T, D>& Vector<T, D>::operator/=(const T& value) {
    #ifdef DEBUG
        if (value == 0) throw std::logic_error("Can't divide by zero");
    #endif

    for (unsigned i = 0; i < D; ++i) _data[i] /= value;
    return *this;
}

template <typename T, unsigned D>
inline Vector<T, D> Vector<T, D>::cross(const Vector<T, D>& value) const {
	static_assert(D >= 3, "Cross product only for dimension 3 (for now anyway..)");

	Vector<T, D> result;
	result(0) = _data[1] * value._data[2] - _data[2] * value._data[1];
	result(1) = _data[2] * value._data[0] - _data[0] * value._data[2];
	result(2) = _data[0] * value._data[1] - _data[1] * value._data[0];
	
	return result;
}

template <typename T, unsigned D>
inline Vector<T, D> Vector<T, D>::unit() const {
    return *this / std::sqrt(this->dot());
}

template <typename T, unsigned D>
inline T Vector<T, D>::dot(const Vector<T, D>& vec) const {
    T result = T(0);
    for (unsigned i = 0; i < D; ++i) result += _data[i] * vec[i];
    return result;
}

template <typename T, unsigned D>
inline T Vector<T, D>::dot() const {
    T result = T(0);
    for (unsigned i = 0; i < D; ++i) result += _data[i] * _data[i];
    return result;
}

template <typename T, unsigned D>
template <typename E>
inline E Vector<T, D>::length() const {
    return std::sqrt(static_cast<E>(dot()));
}


// Definition of non-memberfunctions: --------------------------------------------
template <typename T, unsigned D>
Vector<T, D> operator+(const Vector<T, D>& vec, const Vector<T, D>& other) {
    Vector<T, D> result;
    for (unsigned i = 0; i < D; ++i) result[i] = vec[i] + other[i];
    return result;
}

template <typename T, unsigned D>
Vector<T, D> operator-(const Vector<T, D>& vec, const Vector<T, D>& other) {
    Vector<T, D> result;
    for (unsigned i = 0; i < D; ++i) result[i] = vec[i] - other[i];
    return result;
}

template <typename T, unsigned D>
Vector<T, D> operator*(const Vector<T, D>& vec, const Vector<T, D>& other) {
    Vector<T, D> result;
    for (unsigned i = 0; i < D; ++i) result[i] = vec[i] * other[i];
    return result;
}

template <typename T, unsigned D>
Vector<T, D> operator/(const Vector<T, D>& vec, const Vector<T, D>& other) {
    #ifdef DEBUG
        for (unsigned i = 0; i < D; ++i)
            if (other._data[i] == 0)
                throw std::logic_error("Can't divide by zero");
    #endif

    Vector<T, D> result;
    for (unsigned i = 0; i < D; ++i) result[i] = vec[i] / other[i];
    return result;
}

template <typename T, unsigned D>
inline Vector<T, D> operator*(const Vector<T, D>& vec, const T& value) {
    Vector<T, D> result;
    for (unsigned i = 0; i < D; ++i) result[i] = vec[i] * value;
    return result;
}

template <typename T, unsigned D>
inline Vector<T, D> operator*(const T& value, const Vector<T, D>& vec) {
    return vec * value;
}

template <typename T, unsigned D>
inline Vector<T, D> operator/(const Vector<T, D>& vec, const T& value) {
    Vector<T, D> result;
    for (unsigned i = 0; i < D; ++i) result[i] = vec[i] / value;
    return result;
}

template <typename T, unsigned D>
inline Vector<T, D> operator/(const T& value, const Vector<T, D>& vec) {
    return vec / value;
}

template <typename T, unsigned D>
std::ostream& operator<<(std::ostream& os, const math::Vector<T, D>& vec) {
    os << "(" << vec.x();
    for (unsigned i = 1; i < D; ++i) os << ", " << vec[i];
    os << ")";

    return os;
}

} // math namespace
