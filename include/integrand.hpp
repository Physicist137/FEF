#pragma once
#include <polynomial.hpp>
#include <math.hpp>
#include <functional>
#include <boost/multiprecision/float128.hpp> 
#include <boost/multiprecision/cpp_bin_float.hpp>


template <typename T>
std::function<T(const T&, const T&)> I_integrand(unsigned i, unsigned j) {
	return [=](const T& r, const T& cos) {
		math::Polynomial<T> i_legendre = math::Polynomial<T>::legendre(i);
		math::Polynomial<T> j_legendre = math::Polynomial<T>::legendre(j);
		int p = static_cast<int>(i+j+2);
		return math::pow(r, -p) * i_legendre.at(cos) * j_legendre.at(cos);
	};
}

template <typename T>
std::function<T(const T&, const T&)> G_integrand(unsigned l) {
	return [=](const T& r, const T& cos) {
		math::Polynomial<T> legendre = math::Polynomial<T>::legendre(l);
		int p = static_cast<int>(l);
		return math::pow(r, -p) * legendre.at(cos) * cos;
	};
}


template <typename T>
std::function<T(const T&)>
An_integrand(const T& radius, const T& height, unsigned n) {
	return [=](const T& x) {
		T ee2 = 1.0 - radius*radius / (height * height);
		T xp = math::pow(x, n);
		return xp*math::sqrt(1 + ee2 * ee2 * x * x * (1 - x*x) / (1 - ee2 * x*x) / (1 - ee2 * x*x));
	};
}

template <typename T>
std::function<T(const T&)>
Gl_spheroidal_integrand(const T& radius, const T& height, unsigned l) {
	if (l % 2 == 0) {
		return [=](const T& x) {return T();};
	}

	return [=](const T& x) {
		unsigned p = (l-1) / 2;

		T pi = 3.1415926535;
		T cte = 2*pi / math::pow(radius, 2*p-1);
		T ee2 = 1.0 - radius*radius / (height * height);
		math::Polynomial<T> l_legendre = math::Polynomial<T>::legendre(2*p+1);
		
		T v1 = cte * x * l_legendre.at(x);
		T v2 = math::pow(1 - ee2 * x * x, static_cast<T>(p) - 0.5);
		T v3 = math::sqrt(1 + ee2 * ee2 * x * x * (1 - x*x) / (1 - ee2 * x*x) / (1 - ee2 * x*x));
		return v1 * v2 * v3;
	};
}

template <typename T>
std::function<T(const T&)>
Iij_spheroidal_integrand(const T& radius, const T& height, unsigned i, unsigned j) {
	return [=](const T& x) {
		T pi = 3.1415926535;
		T cte = 2*pi / math::pow(radius, i+j);
		T ee2 = 1.0 - radius*radius / (height * height);
		math::Polynomial<T> i_legendre = math::Polynomial<T>::legendre(i);
		math::Polynomial<T> j_legendre = math::Polynomial<T>::legendre(j);
		
		T v1 = cte * i_legendre.at(x) * j_legendre.at(x);
		T v2 = math::pow(1 - ee2 * x * x, static_cast<int>((i+j) / 2));
		T v3 = math::sqrt(1 + ee2 * ee2 * x * x * (1 - x*x) / (1 - ee2 * x*x) / (1 - ee2 * x*x));
		return v1 * v2 * v3;
	};
}

template <typename T>
std::function<T(const T&)>
Gl_spheroidal_integrand_theta(
const T& radius, const T& height, const std::function<T(const T&)>& r, unsigned l) {
	return [=](const T& theta) {
		T pi = 3.1415926535;
		T cte = 2*pi;
		T ee2 = 1.0 - radius * radius / (height * height);
		T eee = math::sqrt(ee2);
		math::Polynomial<T> l_legendre = math::Polynomial<T>::legendre(l);

		T rv = r(theta);
		T cos = math::cos(theta);
		T sin = math::sin(theta);
		T ratio = rv / radius * eee;

		T sq = math::sqrt(1.0 + ratio * ratio * ratio * ratio * cos * cos * sin * sin);
		T po = math::pow(rv, static_cast<int>(-l+2));
		T ot = cos * sin * l_legendre.at(cos);

		return cte * po * ot * sq;
	};
}

template <typename T>
std::function<T(const T&)>
Iij_spheroidal_integrand_theta(
const T& radius, const T& height, const std::function<T(const T&)>& r, unsigned i, unsigned j) {
	return [=](const T& theta) {
		T pi = 3.1415926535;
		T cte = 2*pi;
		T ee2 = 1.0 - radius * radius / (height * height);
		T eee = math::sqrt(ee2);
		math::Polynomial<T> i_legendre = math::Polynomial<T>::legendre(i);
		math::Polynomial<T> j_legendre = math::Polynomial<T>::legendre(j);

		T rv = r(theta);
		T cos = math::cos(theta);
		T sin = math::sin(theta);
		T ratio = rv / radius * eee;

		T sq = math::sqrt(1.0 + ratio * ratio * ratio * ratio * cos * cos * sin * sin);
		T po = math::pow(rv, static_cast<int>(-i-j));
		T ot = sin * i_legendre.at(cos) * j_legendre.at(cos);

		return cte * po * ot * sq;
	};
}


template <typename T>
std::function<T(const T&)>
Gl_suspended_hemispherical_integrand(const T& radius, const T& length, unsigned l) {
	return [=](const T& x) {
		int p = static_cast<int>(l);

		T pi = 3.1415926535;

		T s;
		if (x > 1) s = 0;
		else s = 1.0 - x*x;

		T sq1 = math::sqrt(radius * radius - length * length * s);

		T r = length * x + sq1;

		math::Polynomial<T> l_legendre = math::Polynomial<T>::legendre(l);

		T factor1 = length * length / r / r * s;
		T factor2 = 1.0 + length * x / sq1;
		T JA = math::sqrt(1.0 + factor1 * factor2 * factor2);
		return 2.0 * pi * math::pow(r, -p+2) * l_legendre.at(x) * x * JA;
	};
}

template <typename T>
std::function<T(const T&)>
Iij_suspended_hemispherical_integrand(const T& radius, const T& length, unsigned i, unsigned j) {
	return [=](const T& x) {
		int p = static_cast<int>(i+j);

		T pi = 3.1415926535;

		T s;
		if (x > 1) s = 0;
		else s = 1.0 - x*x;

		T sq1 = math::sqrt(radius * radius - length * length * s);

		T r = length * x + sq1;

		math::Polynomial<T> i_legendre = math::Polynomial<T>::legendre(i);
		math::Polynomial<T> j_legendre = math::Polynomial<T>::legendre(j);

		T factor1 = length * length / r / r * s;
		T factor2 = 1.0 + length * x / sq1;
		T JA = math::sqrt(1.0 + factor1 * factor2 * factor2);
		return 2.0 * pi * math::pow(r, -p) * i_legendre.at(x) * j_legendre.at(x) * JA;
	};
}

template <typename T>
std::function<T(const T&)>
Gl_cylindrical_integrand(const T& radius, const T& length, unsigned l) {
	return [=](const T& x) {
		int p = static_cast<int>(l);

		T pi = 3.1415926535;

		T s;
		if (x > 1) s = 0;
		else s = 1.0 - x*x;

		T sqr = math::sqrt(1.0 - x*x);
		T area = math::pow(sqr, p-3);
		math::Polynomial<T> l_legendre = math::Polynomial<T>::legendre(l);
		return 2.0 * pi * math::pow(radius, -p+2) * l_legendre.at(x) * x * area;
	};
}

template <typename T>
std::function<T(const T&)>
Iij_cylindrical_integrand(const T& radius, const T& length, unsigned i, unsigned j) {
	return [=](const T& x) {
		int p = static_cast<int>(i+j);

		T pi = 3.1415926535;

		T s;
		if (x > 1) s = 0;
		else s = 1.0 - x*x;

		T sqr = math::sqrt(1.0 - x*x);
		T area = math::pow(sqr, p-1);
		math::Polynomial<T> i_legendre = math::Polynomial<T>::legendre(i);
		math::Polynomial<T> j_legendre = math::Polynomial<T>::legendre(j);
		return 2.0 * pi * math::pow(radius, -p) * i_legendre.at(x) * j_legendre.at(x) * area;
	};
}

