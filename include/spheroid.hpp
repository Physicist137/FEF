#pragma once
#include <integrand.hpp>
#include <polynomial.hpp>
#include <cte.hpp>

template <typename T>
T Gl_spheroidal_summation_approx(float r, float h, unsigned l) {
	if (l % 2 == 0) return T(0);

	unsigned p = (l-1) / 2;
	T radius = T(r);
	T height = T(h);

	T pt = static_cast<T>(p);
	T pi = T(math::cte::pi_string);
	T cte = 2 * pi / boost::multiprecision::pow(radius, 2*p-1);
	T ee2 = 1.0 - radius * radius / (height * height);
	math::Polynomial<T> p_legendre = math::Polynomial<T>::legendre(2*p+1);
	T result = 0.0;
	T bin = 1.0;

	for (unsigned n = 0; n < 10; ++n) {
		T mul = 1;
		T nt = static_cast<T>(n);

		if (n >= 1) bin *= (pt - 0.5 - nt+1.0) / nt;
		if (n % 2 == 0) mul = 1;
		else mul = -1; 

		for (unsigned k = 0; k <= p; ++k) {
			T kt = static_cast<T>(k);
			T main = p_legendre[2*k+1] * bin * boost::multiprecision::pow(ee2, n);
			T moment = 1.0 / (nt+kt+1.0+0.5) + ee2 * ee2 / 2.0 / (nt+kt+1.0+1.5);

			std::cout << "(" << pt-0.5 << ";" << nt << ") = " << bin << std::endl;
			

			result += mul * main * moment;
		}

	}

	return cte * result;
}

template <typename T>
T Gl_spheroidal_summation_exact(float r, float h, unsigned l, const std::vector<T>& An = {}) {
	if (l % 2 == 0) return T(0);

	unsigned p = (l-1) / 2;
	T radius = T(r);
	T height = T(h);

	T pt = static_cast<T>(p);
	T pi = T(math::cte::pi_string);
	T cte = 2 * pi / boost::multiprecision::pow(radius, 2*p-1);
	T ee2 = 1.0 - radius * radius / (height * height);
	math::Polynomial<T> p_legendre = math::Polynomial<T>::legendre(2*p+1);
	T result = 0.0;
	T partial_sum = 0.0;
	T bin = 1.0;
	T error = 0.0;

	unsigned n = 0;
	do {
		T mul = 1;
		T nt = static_cast<T>(n);

		if (n % 2 == 0) mul = 1;
		else mul = -1; 

		for (unsigned k = 0; k <= p; ++k) {
			T main = p_legendre[2*k+1] * bin;
			T pe = boost::multiprecision::pow(ee2, n);

			T moment;
			if (An.empty()) moment = integrate<T>(-1.0, 1.0, An_integrand<T>(radius, height, 2*(k+n+1)));
			else if (An.size() <= 2*(k+n+1)) moment = integrate<T>(-1.0, 1.0, An_integrand<T>(radius, height, 2*(k+n+1)));
			else moment = An[2*(k+n+1)];

			partial_sum += mul * main * (moment * pe);
		}

		bin *= (pt - 0.5 - nt) / (nt + 1.0);
		result += partial_sum;
		error = boost::multiprecision::abs(partial_sum / result);
		partial_sum = 0.0;
		++n;
	} while (2*p > n+1  or  error > 1e-5);

	return cte * result;
}

template <typename T>
T Iij_spheroidal_summation_approx(float r, float h, unsigned i, unsigned j) {
	if ((i+j) % 2 == 1) return T(0);

	T radius = static_cast<T>(r);
	T height = static_cast<T>(h);

	T ij = static_cast<T>((i+j) / 2);
	T pi = T(math::cte::pi_string);
	T cte = 2.0 * pi / boost::multiprecision::pow(radius, i+j);
	T ee2 = 1.0 - radius * radius / (height * height);
	math::Polynomial<T> i_legendre = math::Polynomial<T>::legendre(i);
	math::Polynomial<T> j_legendre = math::Polynomial<T>::legendre(j);
	T result = 0.0;
	T bin = 1.0;

	for (unsigned n = 0; n <= (i+j) / 2; ++n) {
		T mul = 1;
		T nt = static_cast<T>(n);

		if (n % 2 == 0) mul = 1;
		else mul = -1; 

		for (unsigned u = 0; u <= i; ++u) {
			for (unsigned v = 0; v <= j; ++v) {
				if ((u + v) % 2 == 1) continue;

				T uv = static_cast<T>((u+v) / 2);
				T main = i_legendre[u] * j_legendre[v] * bin * boost::multiprecision::pow(ee2, n);
				T moment = 1.0 / (nt+uv+0.5) + ee2 * ee2 / 2.0 / (nt+uv+1.5);

				std::cout << "(" << ij << ";" << nt << ") = " << bin << std::endl;

				result += mul * main * moment;
			}
		}

		bin *= (ij - nt) / (nt + 1.0);
	}

	return cte * result;
}

template <typename T>
T Iij_spheroidal_summation_exact(float r, float h, unsigned i, unsigned j, const std::vector<T>& An = {}) {
	if ((i+j) % 2 == 1) return T(0);

	T radius = static_cast<T>(r);
	T height = static_cast<T>(h);

	T ij = static_cast<T>((i+j) / 2);
	T pi = T(math::cte::pi_string);
	T cte = 2 * pi / boost::multiprecision::pow(radius, i+j);
	T ee2 = 1.0 - radius * radius / (height * height);
	math::Polynomial<T> i_legendre = math::Polynomial<T>::legendre(i);
	math::Polynomial<T> j_legendre = math::Polynomial<T>::legendre(j);
	T result = 0.0;
	T bin = 1.0;

	for (unsigned n = 0; n <= (i+j) / 2; ++n) {
		T mul = 1;
		T nt = static_cast<T>(n);

		if (n >= 1) bin *= (ij - nt+1.0) / nt;
		if (n % 2 == 0) mul = 1;
		else mul = -1; 

		for (unsigned u = 0; u <= i; ++u) {
			for (unsigned v = 0; v <= j; ++v) {
				T main = i_legendre[u] * j_legendre[v] * bin * boost::multiprecision::pow(ee2, n);

				T moment;
				if (An.empty()) moment = integrate<T>(-1, 1, An_integrand<T>(radius, height, u+v+2*n));
				else if (An.size() <= u+v+2*n) moment = integrate<T>(-1, 1, An_integrand<T>(radius, height, u+v+2*n));
				else moment = An[u+v+2*n];

				// std::cout << "(" << ij << ";" << nt << ") = " << bin << std::endl;

				result += mul * main * moment;
			}
		}
	}

	return cte * result;
}

