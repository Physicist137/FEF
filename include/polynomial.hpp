#pragma once
#include <vector>
#include <cmath>
#include <initializer_list>
#include <ostream>
#include <algorithm>

/*!
This is a polynomial over a field, that is, the ring Field[x].
The coeficients starts with the most non important term:
_coeff[0] * x^0 + _coeff[1] * x^1 + ... + _coeff[n] * x^n
*/

/*
OTHER FEATURES TO BE DONE:
	- To calculate roots of the polynomials.
	- More static functions for more polynomials (not noly legendre)
*/

namespace math {
template <typename Field>
class Polynomial {
	std::vector<Field> _coeff;

public:
    Polynomial() : _coeff{Field()} {}
	Polynomial(const Field& value) : _coeff{value} {}
    Polynomial(std::initializer_list<Field> list) : _coeff(list) {}

    Polynomial(const Polynomial<Field>& other) : _coeff(other._coeff) {}
    Polynomial(Polynomial<Field>& other) : _coeff(other._coeff) {}

	inline unsigned degree() const {return _coeff.size() - 1;}
	inline unsigned size() const {return _coeff.size();}
	
	Polynomial& push_back(const Field& value);

    const Field& operator[](unsigned i) const;
    Field& operator[](unsigned i);

    const Field& operator()(unsigned i) const;
    Field& operator()(unsigned i);

    Polynomial operator+() const;
    Polynomial operator-() const;
    Polynomial operator/(const Field& value);

    Polynomial& operator+=(const Polynomial& other);
    Polynomial& operator-=(const Polynomial& other);
    Polynomial& operator*=(const Polynomial& other);

    Polynomial& operator+=(const Field& other);
    Polynomial& operator-=(const Field& other);
    Polynomial& operator*=(const Field& value);
    Polynomial& operator/=(const Field& value);

	Polynomial& operator=(const Polynomial& other);
	bool operator==(const Polynomial& other) const;
	bool operator!=(const Polynomial& other) const;

	/// Evaluates the polynomial at x=value.
	Field evaluateAt(const Field& value) const;
	inline Field at(const Field& value) const {return evaluateAt(value);}

public:
	/// 1 and x polynomials.
	static Polynomial one();
	static Polynomial x();

	/// Generate Legendre polynomials of some given order.
	static Polynomial legendre(unsigned order);
	static std::vector<Polynomial> legendre_base(unsigned order);
};

// Declaration of functions: -----------------------------------------------------
template <typename Field> Polynomial<Field> operator+(const Polynomial<Field>& vec, const Polynomial<Field>& other);
template <typename Field> Polynomial<Field> operator-(const Polynomial<Field>& vec, const Polynomial<Field>& other);
template <typename Field> Polynomial<Field> operator*(const Polynomial<Field>& vec, const Polynomial<Field>& other);

template <typename Field> Polynomial<Field> operator+(const Polynomial<Field>& poly, const Field& value);
template <typename Field> Polynomial<Field> operator-(const Polynomial<Field>& poly, const Field& value);
template <typename Field> Polynomial<Field> operator+(const Field& value, const Polynomial<Field>& poly);
template <typename Field> Polynomial<Field> operator-(const Field& value, const Polynomial<Field>& poly);

template <typename Field> Polynomial<Field> operator*(const Polynomial<Field>& vec, const Field& value);
template <typename Field> Polynomial<Field> operator*(const Field& value, const Polynomial<Field>& vec);


// Definition of memberfunctions: ------------------------------------------------
template <typename Field>
Polynomial<Field>& Polynomial<Field>::push_back(const Field& value) {
	_coeff.push_back(value);
	return *this;
}

template <typename Field>
inline const Field& Polynomial<Field>::operator[](unsigned i) const {
	return _coeff[i];
}

template <typename Field>
inline Field& Polynomial<Field>::operator[](unsigned i) {
	return _coeff[i];
}

template <typename Field>
inline const Field& Polynomial<Field>::operator()(unsigned i) const {
	return this->operator[](i);
}

template <typename Field>
inline Field& Polynomial<Field>::operator()(unsigned i) {
	return this->operator[]();
}

template <typename Field>
inline Polynomial<Field> Polynomial<Field>::operator+() const {
	return *this;
}

template <typename Field>
inline Polynomial<Field> Polynomial<Field>::operator-() const {
	Polynomial<Field> result;
	unsigned size = _coeff.size();

	result[0] = -_coeff[0];
	for (unsigned i = 1; i < size; ++i) result._coeff.push_back(-_coeff[i]);
	return result;
}

template <typename Field>
inline Polynomial<Field> Polynomial<Field>::operator/(const Field& value) {
	Polynomial<Field> result;
	unsigned size = _coeff.size();
	
	result[0] = _coeff[0] / value;
	for (unsigned i = 1; i < size; ++i) result._coeff.push_back(_coeff[i] / value);
	return result;
}

template <typename Field>
Polynomial<Field>& Polynomial<Field>::operator+=(const Polynomial<Field>& other) {
	unsigned size_other = other._coeff.size();
	unsigned size_this = _coeff.size();
	unsigned size_min = std::min(size_other, size_this);

	// Make the sum from the existing coeficients.
	for (unsigned i = 0; i < size_min; ++i)	_coeff[i] += other[i];

	// Create more coeficients if necessary.
	if (size_other > size_this) {
		for (unsigned i = size_min; i < size_other; ++i) _coeff.push_back(other[i]);
	}

	return *this;
}

template <typename Field>
Polynomial<Field>& Polynomial<Field>::operator-=(const Polynomial<Field>& other) {
	unsigned size_other = other._coeff.size();
	unsigned size_this = _coeff.size();
	unsigned size_min = std::min(size_other, size_this);

	// Make the sum from the existing coeficients.
	for (unsigned i = 0; i < size_min; ++i)	_coeff[i] -= other[i];

	// Create more coeficients if necessary.
	if (size_other > size_this) {
		for (unsigned i = size_min; i < size_other; ++i) _coeff.push_back(-other[i]);
	}

	return *this;
}

template <typename Field>
Polynomial<Field>& Polynomial<Field>::operator*=(const Polynomial<Field>& other) {
	return *this = (*this) * other;
}


template <typename Field>
Polynomial<Field>& Polynomial<Field>::operator+=(const Field& value) {
	_coeff[0] += value;
	return *this;
}

template <typename Field>
Polynomial<Field>& Polynomial<Field>::operator-=(const Field& value) {
	_coeff[0] -= value;
	return *this;
}

template <typename Field>
Polynomial<Field>& Polynomial<Field>::operator*=(const Field& value) {
	unsigned size = _coeff.size();
	for (unsigned i = 0; i < size; ++i) _coeff[i] *= value;
	return *this;
}

template <typename Field>
Polynomial<Field>& Polynomial<Field>::operator/=(const Field& value) {
	unsigned size = _coeff.size();
	for (unsigned i = 0; i < size; ++i) _coeff[i] /= value;
	return *this;
}

template <typename Field>
Polynomial<Field>& Polynomial<Field>::operator=(const Polynomial& other) {
	unsigned size = other.size();
	
	_coeff.clear();
	for (unsigned i = 0; i < size; ++i) _coeff.push_back(other[i]);
	return *this;
}

template <typename Field>
bool Polynomial<Field>::operator==(const Polynomial& other) const {
	unsigned size = _coeff.size();
	if (size != other.size()) return false;
	
	for (unsigned i = 0; i < size; ++i)
		if (_coeff[i] != other[i]) return false;

	return true;
}

template <typename Field>
bool Polynomial<Field>::operator!=(const Polynomial& other) const {
	return not (*this == other);
}


template <typename Field>
Field Polynomial<Field>::evaluateAt(const Field& value) const {
	unsigned size = _coeff.size();

	// Compute the x^0
	Field xpower = 1.0;
	Field result = _coeff[0];

	// Compute the x^n part.
	for (unsigned n = 1; n < size; ++n) {
		xpower *= value;
		result += _coeff[n] * xpower;
	}

	return result;
}


// Static Member functions -----------------------------------------------------------

template <typename Field>
inline Polynomial<Field> Polynomial<Field>::one() {
	return Polynomial<Field>{
		static_cast<Field>(1.0)
	};
}

template <typename Field>
inline Polynomial<Field> Polynomial<Field>::x() {
	return Polynomial<Field>{
		static_cast<Field>(0.0),
		static_cast<Field>(1.0)
	};
}

template <typename Field>
Polynomial<Field> Polynomial<Field>::legendre(unsigned order) {
	Polynomial<Field> one = Polynomial<Field>::one();
	Polynomial<Field> x = Polynomial<Field>::x();

	if (order == 0) return one;
	if (order == 1) return x;

	Polynomial<Field> previous = one;
	Polynomial<Field> current = x;
	for (unsigned n = 1; n < order; ++n) {
		Field fn = static_cast<Field>(n);
		Polynomial<Field> next = ((2*fn+1) * x * current - fn * previous) / (fn + 1);
		previous = current;
		current = next;
	}

	return current;
}

template <typename Field>
std::vector<Polynomial<Field>> Polynomial<Field>::legendre_base(unsigned order) {
	std::vector<Polynomial<Field>> result;
	Polynomial<Field> x = {0.0, 1.0};

	result.push_back(Polynomial<Field>{1.0});
	if (order == 0) return result;

	result.push_back(x);
	if (order == 1) return result;

	Polynomial<Field> previous = 1.0;
	Polynomial<Field> current = x;
	for (unsigned n = 1; n < order; ++n) {
		Field fn = static_cast<Field>(n);
		
		Polynomial<Field> next = ((2*fn+1) * x * current - fn * previous) / (fn + 1);
		result.push_back(next);

		previous = current;
		current = next;
	}

	return result;
}


// -----------------------------------------------------------------------------------
template <typename Field>
Polynomial<Field> operator+(const Polynomial<Field>& poly, const Polynomial<Field>& other) {
	Polynomial<Field> result;
	unsigned poly_size = poly.size();
	unsigned other_size = other.size();
	unsigned min_size = std::min(poly_size, other_size);

	result[0] = poly[0] + other[0];
	for (unsigned i = 1; i < min_size; ++i) result.push_back(poly[i] + other[i]);

	if (poly_size > other_size)
		for (unsigned i = min_size; i < poly_size; ++i) result.push_back(poly[i]);

	if (other_size > poly_size)
		for (unsigned i = min_size; i < other_size; ++i) result.push_back(other[i]);

	return result;
}

template <typename Field>
Polynomial<Field> operator-(const Polynomial<Field>& poly, const Polynomial<Field>& other) {
	return poly + (-other);
}

template <typename Field>
Polynomial<Field> operator*(const Polynomial<Field>& poly, const Polynomial<Field>& other) {
	Polynomial<Field> result;
	unsigned poly_size = poly.size();
	unsigned other_size = other.size();
	unsigned result_size = poly.degree() + other.degree() + 1;

	for (unsigned n = 1; n < result_size; ++n)
		result.push_back(0.0);

	for (unsigned n = 0; n < poly_size; ++n) {
		for (unsigned k = 0; k < other_size; ++k) {
			result[n+k] += poly[n] * other[k];
		}
	}

	return result;
}


template <typename Field>
Polynomial<Field> operator+(const Polynomial<Field>& poly, const Field& value) {
	Polynomial<Field> result(poly);
	result[0] += value;
	return result;
}

template <typename Field>
Polynomial<Field> operator-(const Polynomial<Field>& poly, const Field& value) {
	return poly + (-value);
}

template <typename Field>
Polynomial<Field> operator+(const Field& value, const Polynomial<Field>& poly) {
	return poly + value;
}

template <typename Field>
Polynomial<Field> operator-(const Field& value, const Polynomial<Field>& poly) {
	return poly - value;
}

template <typename Field>
Polynomial<Field> operator*(const Polynomial<Field>& poly, const Field& value) {
	Polynomial<Field> result;
	unsigned size = poly.size();

	result[0] = poly[0] * value;
	for (unsigned i = 1; i < size; ++i) result.push_back(poly[i] * value);

	return result;
}

template <typename Field>
Polynomial<Field> operator*(const Field& value, const Polynomial<Field>& poly) {
	return poly * value;
}


} // Namespace math.


template <typename Field>
std::ostream& operator<<(std::ostream& os, const math::Polynomial<Field>& polynomial) {
	unsigned size = polynomial.size();

	os << "[" << polynomial[0];
	for (unsigned i = 1; i < size; ++i) os << ", " << polynomial[i];
	os << "]";
	
	return os;
}
